#!/usr/bin/env python

import logging
import json
import sys
import time
from typing import Union
from multiprocessing import cpu_count

from Bio.Align import MultipleSeqAlignment
from .args_processing import process_args
from .exceptions import InvalidInputFileFormat
from .files import (
    get_alignment_and_format,
    FileFormat,
    write_debug_log_file,
    get_custom_sites_to_trim,
)
from .helpers import (
    create_msa,
    get_seq_type,
    get_gap_chars,
    write_msa,
    write_complement,
    SeqType,
)
from .logger import logger, log_file_logger
from .modes import TrimmingMode
from .msa import MSA
from .parser import create_parser
from .settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS
from .smart_gap_helper import smart_gap_threshold_determination
from .version import __version__ as current_version
from .warnings import (
    warn_if_all_sites_were_trimmed,
    warn_if_entry_contains_only_gaps,
)
from .write import (
    write_user_args,
    write_output_stats,
    write_output_files_message,
)

from dataclasses import dataclass


@dataclass
class TrimRun:
    alignment: MultipleSeqAlignment
    msa: MSA
    gap_characters: list
    sequence_type: SeqType
    input_file_format: FileFormat
    output_file_format: FileFormat
    gaps: float
    codon: bool
    version: str = current_version

    @property
    def complement(self):
        return self.msa.complement_to_bio_msa()

    @property
    def trimmed(self):
        return self.msa.to_bio_msa()


KPI_FAMILY_MODES = {
    TrimmingMode.kpi,
    TrimmingMode.kpi_gappy,
    TrimmingMode.kpi_smart_gap,
    TrimmingMode.kpic,
    TrimmingMode.kpic_gappy,
    TrimmingMode.kpic_smart_gap,
}


def determine_effective_threads(
    requested_threads: int,
    mode: TrimmingMode,
    n_sequences: int,
    alignment_length: int,
) -> int:
    """
    Heuristic thread selection.

    KPI/KPIC-family workloads can regress with high thread counts due to task
    scheduling overhead and relatively small per-task compute.
    """
    if requested_threads <= 1:
        return 1

    max_threads = max(1, cpu_count() or 1)
    requested_threads = min(requested_threads, max_threads)
    workload = n_sequences * alignment_length

    if mode in KPI_FAMILY_MODES:
        # Keep single-threaded for smaller workloads where parallel overhead
        # commonly dominates.
        if workload < 100_000_000:
            return 1
        # Moderate workloads usually benefit from modest concurrency.
        if workload < 250_000_000:
            return min(requested_threads, 2)
        # Larger workloads can use a few more workers.
        if workload < 500_000_000:
            return min(requested_threads, 4)

    return requested_threads


def run(
    input_file: str,
    input_file_format: FileFormat,
    output_file: str,
    output_file_format: FileFormat,
    auxiliary_file: str,
    sequence_type: Union[SeqType, None],
    gaps: float,
    gap_characters: Union[list, None],
    complement: bool,
    codon: bool,
    mode: TrimmingMode,
    use_log: bool,
    quiet: bool,
    ends_only: bool,
    threads: int = 1,
):
    alignment, input_file_format = get_alignment_and_format(input_file, input_file_format)

    if threads < 1:
        raise ValueError("threads must be an integer >= 1")

    sequence_type = sequence_type or get_seq_type(alignment)

    if not gap_characters:
        gap_characters = get_gap_chars(sequence_type)

    if not output_file_format:
        output_file_format = input_file_format
    else:
        output_file_format = FileFormat(output_file_format)

    site_positions_to_trim = None
    if mode == TrimmingMode.cst:
        aln_length = alignment.get_alignment_length()
        site_positions_to_trim = (
            get_custom_sites_to_trim(auxiliary_file, aln_length) or []
        )

    effective_threads = determine_effective_threads(
        threads,
        mode,
        len(alignment),
        alignment.get_alignment_length(),
    )

    msa = create_msa(alignment, gap_characters, effective_threads)

    # determine smart_gap threshold
    if mode in {
        TrimmingMode.smart_gap,
        TrimmingMode.kpi_smart_gap,
        TrimmingMode.kpic_smart_gap,
    }:
        gaps = smart_gap_threshold_determination(
            alignment,
            gap_characters,
            seq_records=msa.seq_records,
        )

    msa.trim(
        mode,
        gap_threshold=gaps,
        site_positions_to_trim=site_positions_to_trim,
        codon=codon,
        ends_only=ends_only,
    )

    trim_run = TrimRun(
        alignment,
        msa,
        gap_characters,
        sequence_type,
        input_file_format,
        output_file_format,
        gaps,
        codon,
    )

    return trim_run, msa.stats


def execute(
    input_file: str,
    input_file_format: FileFormat,
    output_file: str,
    output_file_format: FileFormat,
    sequence_type: Union[SeqType, None],
    gaps: float,
    gap_characters: Union[list, None],
    complement: bool,
    codon: bool,
    ends_only: bool,
    mode: TrimmingMode,
    use_log: bool,
    quiet: bool,
    dry_run: bool = False,
    validate_only: bool = False,
    report_json: Union[str, None] = None,
    auxiliary_file: str = None,
    threads: int = 1,
    **kwargs,
) -> None:
    fh = None
    original_logger_disabled = logger.disabled
    original_log_level = log_file_logger.level
    original_log_propagate = log_file_logger.propagate
    try:
        if use_log and not dry_run and not validate_only:
            log_file_logger.setLevel(logging.DEBUG)
            log_file_logger.propagate = False
            fh = logging.FileHandler(f"{output_file}.log", mode="w")
            fh.setLevel(logging.DEBUG)
            log_file_logger.addHandler(fh)

        if quiet:
            logger.disabled = True

        start_time = time.time()

        if validate_only:
            try:
                alignment, detected_input_format = get_alignment_and_format(
                    input_file, input_file_format
                )
            except InvalidInputFileFormat:
                logger.error(
                    f"""Format type could not be read.\nPlease check acceptable input file formats: {", ".join([file_format.value for file_format in FileFormat])}"""
                )
                return

            validated_sequence_type = sequence_type or get_seq_type(alignment)
            validated_gap_characters = gap_characters or get_gap_chars(
                validated_sequence_type
            )
            resolved_output_format = (
                detected_input_format
                if not output_file_format
                else FileFormat(output_file_format)
            )

            if mode == TrimmingMode.cst:
                aln_length = alignment.get_alignment_length()
                get_custom_sites_to_trim(auxiliary_file, aln_length)

            logger.info("Validation successful. No trimming was performed.")

            if report_json:
                _write_report_json(
                    report_json,
                    dict(
                        status="validated",
                        validate_only=True,
                        dry_run=False,
                        input_file=input_file,
                        output_file=output_file,
                        input_file_format=detected_input_format.value,
                        output_file_format=resolved_output_format.value,
                        sequence_type=validated_sequence_type.value,
                        gaps=gaps,
                        gap_characters=validated_gap_characters,
                        mode=mode.value,
                        codon=codon,
                        complement=complement,
                        ends_only=ends_only,
                        threads_requested=threads,
                        runtime_seconds=round(time.time() - start_time, 6),
                    ),
                )
            return

        try:
            trim_run, stats = run(
                input_file,
                input_file_format,
                output_file,
                output_file_format,
                auxiliary_file,
                sequence_type,
                gaps,
                gap_characters,
                complement,
                codon,
                mode,
                use_log,
                quiet,
                ends_only,
                threads,
            )
        except InvalidInputFileFormat:
            logger.error(
                f"""Format type could not be read.\nPlease check acceptable input file formats: {", ".join([file_format.value for file_format in FileFormat])}"""
            )
            return

        # display to user what args are being used in stdout
        write_user_args(
            input_file,
            trim_run.input_file_format,
            output_file,
            trim_run.output_file_format,
            trim_run.sequence_type,
            trim_run.gaps,
            trim_run.gap_characters,
            mode,
            complement,
            codon,
            use_log,
            ends_only,
        )

        if dry_run:
            logger.info("Dry run enabled. No output files were written.")
        else:
            write_output_files_message(output_file, complement, use_log)

        if (not dry_run) and use_log:
            warn_if_all_sites_were_trimmed(trim_run.msa)
            warn_if_entry_contains_only_gaps(trim_run.msa)
            write_debug_log_file(trim_run.msa)

        base_metadata = trim_run.alignment.annotations.get("ecomp_metadata")

        if not dry_run:
            write_msa(
                trim_run.msa,
                output_file,
                trim_run.output_file_format,
                base_metadata=base_metadata,
            )

            # if the -c/--complementary argument was used, create an alignment of the trimmed sequences
            if complement:
                write_complement(
                    trim_run.msa,
                    output_file,
                    trim_run.output_file_format,
                    base_metadata=base_metadata,
                )

        if report_json:
            _write_report_json(
                report_json,
                dict(
                    status="completed",
                    validate_only=False,
                    dry_run=dry_run,
                    input_file=input_file,
                    output_file=output_file,
                    input_file_format=trim_run.input_file_format.value,
                    output_file_format=trim_run.output_file_format.value,
                    sequence_type=trim_run.sequence_type.value,
                    gaps=trim_run.gaps,
                    gap_characters=trim_run.gap_characters,
                    mode=mode.value,
                    codon=codon,
                    complement=complement,
                    ends_only=ends_only,
                    threads_requested=threads,
                    threads_effective=trim_run.msa._threads,
                    stats=stats.summary,
                    runtime_seconds=round(time.time() - start_time, 6),
                ),
            )

        write_output_stats(stats, start_time)
    finally:
        if fh is not None:
            log_file_logger.removeHandler(fh)
            fh.close()
        logger.disabled = original_logger_disabled
        log_file_logger.setLevel(original_log_level)
        log_file_logger.propagate = original_log_propagate


def _write_report_json(path: str, payload: dict) -> None:
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def main(argv=None):
    """
    Function that parses and collects arguments
    """
    parser = create_parser()
    argv = sys.argv[1:] if argv is None else argv

    if len(argv) == 0:
        parser.print_help(sys.stderr)
        return

    args = parser.parse_args(argv)

    execute(**process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])

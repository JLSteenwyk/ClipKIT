#!/usr/bin/env python

import logging
import sys
import time
from typing import Union

from Bio.Align import MultipleSeqAlignment
from .args_processing import process_args
from .exceptions import InvalidInputFileFormat
from .files import get_alignment_and_format, FileFormat
from .helpers import (
    get_seq_type,
    get_gap_chars,
    keep_trim_and_log,
    write_keep_msa,
    write_trim_msa,
    SeqType,
)
from .logger import logger, log_file_logger
from .modes import TrimmingMode
from .msa import MSA
from .parser import create_parser
from .settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS
from .smart_gap_helper import smart_gap_threshold_determination
from .stats import TrimmingStats
from .warnings import (
    warn_if_all_sites_were_trimmed,
    warn_if_entry_contains_only_gaps,
)
from .write import (
    write_determining_smart_gap_threshold,
    write_user_args,
    write_output_stats,
)

from dataclasses import dataclass


@dataclass
class TrimRun:
    alignment: MultipleSeqAlignment
    keep_msa: MSA
    trim_msa: MSA
    gap_characters: list
    sequence_type: SeqType
    input_file_format: FileFormat
    output_file_format: FileFormat
    site_classification_counts: dict
    gaps: float


def run(
    input_file: str,
    input_file_format: FileFormat,
    output_file: str,
    output_file_format: FileFormat,
    sequence_type: Union[SeqType, None],
    gaps: float,
    gap_characters: Union[list, None],
    complement: bool,
    mode: TrimmingMode,
    use_log: bool,
    quiet: bool,
):
    try:
        alignment, input_file_format = get_alignment_and_format(
            input_file, input_file_format
        )
    except InvalidInputFileFormat:
        return logger.error(
            f"""Format type could not be read.\nPlease check acceptable input file formats: {", ".join([file_format.value for file_format in FileFormat])}"""
        )

    sequence_type = sequence_type or get_seq_type(alignment)

    if not gap_characters:
        gap_characters = get_gap_chars(sequence_type)

    if not output_file_format:
        output_file_format = input_file_format
    else:
        output_file_format = FileFormat(output_file_format)

    # determine smart_gap threshold
    if mode in {
        TrimmingMode.smart_gap,
        TrimmingMode.kpi_smart_gap,
        TrimmingMode.kpic_smart_gap,
    }:
        write_determining_smart_gap_threshold()
        gaps = smart_gap_threshold_determination(alignment, gap_characters, quiet)

    # instantiates MSAs to track what we keep/trim from the alignment
    keep_msa, trim_msa, site_classification_counts = keep_trim_and_log(
        alignment, gaps, mode, use_log, output_file, complement, gap_characters, quiet
    )

    return TrimRun(
        alignment,
        keep_msa,
        trim_msa,
        gap_characters,
        sequence_type,
        input_file_format,
        output_file_format,
        site_classification_counts,
        gaps,
    )


def execute(
    input_file: str,
    input_file_format: FileFormat,
    output_file: str,
    output_file_format: FileFormat,
    sequence_type: Union[SeqType, None],
    gaps: float,
    gap_characters: Union[list, None],
    complement: bool,
    mode: TrimmingMode,
    use_log: bool,
    quiet: bool,
    **kwargs,
) -> None:
    if use_log:
        log_file_logger.setLevel(logging.DEBUG)
        log_file_logger.propagate = False
        fh = logging.FileHandler(f"{output_file}.log", mode="w")
        fh.setLevel(logging.DEBUG)
        log_file_logger.addHandler(fh)

    if quiet:
        logger.disabled = True

    # for reporting runtime duration to user
    start_time = time.time()

    trim_run = run(
        input_file,
        input_file_format,
        output_file,
        output_file_format,
        sequence_type,
        gaps,
        gap_characters,
        complement,
        mode,
        use_log,
        quiet,
    )

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
        use_log,
    )

    if use_log:
        warn_if_all_sites_were_trimmed(trim_run.keep_msa)
        warn_if_entry_contains_only_gaps(trim_run.keep_msa, trim_run.sequence_type)

    write_keep_msa(trim_run.keep_msa, output_file, trim_run.output_file_format)

    # if the -c/--complementary argument was used, create an alignment of the trimmed sequences
    if complement:
        write_trim_msa(trim_run.trim_msa, output_file, trim_run.output_file_format)

    stats = TrimmingStats(trim_run.alignment, trim_run.keep_msa, trim_run.trim_msa)
    write_output_stats(stats, start_time)


def main(argv=None):
    """
    Function that parses and collects arguments
    """
    parser = create_parser()
    args = parser.parse_args()

    execute(**process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])

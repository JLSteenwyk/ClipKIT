#!/usr/bin/env python

import logging
import sys
import time
from typing import Union

from Bio.Align import MultipleSeqAlignment
from .args_processing import process_args
from .exceptions import InvalidInputFileFormat
from .files import get_alignment_and_format, FileFormat, write_debug_log_file
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
    write_processing_aln,
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


def run(
    input_file: str,
    input_file_format: FileFormat,
    output_file: str,
    output_file_format: FileFormat,
    sequence_type: Union[SeqType, None],
    gaps: float,
    gap_characters: Union[list, None],
    complement: bool,
    codon: bool,
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
        gaps = smart_gap_threshold_determination(alignment, gap_characters)

    msa = create_msa(alignment, gap_characters)
    msa.trim(mode, gap_threshold=gaps, site_positions_to_trim=None, codon=codon)

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

    write_processing_aln()
    trim_run, stats = run(
        input_file,
        input_file_format,
        output_file,
        output_file_format,
        sequence_type,
        gaps,
        gap_characters,
        complement,
        codon,
        mode,
        use_log,
        quiet,
    )
    write_output_files_message(output_file, complement, use_log)

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
    )

    if use_log:
        warn_if_all_sites_were_trimmed(trim_run.msa)
        warn_if_entry_contains_only_gaps(trim_run.msa)
        write_debug_log_file(trim_run.msa)

    write_msa(trim_run.msa, output_file, trim_run.output_file_format)

    # if the -c/--complementary argument was used, create an alignment of the trimmed sequences
    if complement:
        write_complement(trim_run.msa, output_file, trim_run.output_file_format)

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

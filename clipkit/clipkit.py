#!/usr/bin/env python

import getopt
import logging
import os.path
import sys
import time
from typing import Union

from .args_processing import process_args
from .exceptions import InvalidInputFileFormat
from .files import get_alignment_and_format, FileFormat
from .helpers import (
    get_seq_type,
    keep_trim_and_log,
    write_keepD,
    write_trimD
)
from .helpers import SeqType
from .logger import logger, log_file_logger
from .modes import TrimmingMode
from .parser import create_parser
from .smart_gap_helper import smart_gap_threshold_determination
from .stats import TrimmingStats
from .warnings import (
    warn_if_all_sites_were_trimmed,
    warn_if_entry_contains_only_gaps,
)
from .write import write_determining_smart_gap_threshold, write_user_args, write_output_stats


def execute(
    input_file: str,
    input_file_format: FileFormat,
    output_file: str,
    output_file_format: FileFormat,
    sequence_type: Union[SeqType, None],
    gaps: float,
    complement: bool,
    mode: TrimmingMode,
    use_log: bool,
    quiet: bool,
):
    if use_log:
        log_file_logger.setLevel(logging.DEBUG)
        log_file_logger.propagate = False
        fh = logging.FileHandler(f"{output_file}.log", mode="w")
        fh.setLevel(logging.DEBUG)
        log_file_logger.addHandler(fh)

    if quiet:
        logger.disabled = True

    # create start time logger
    start_time = time.time()

    # read in alignment and save the format of the alignment
    try:
        alignment, input_file_format = get_alignment_and_format(
            input_file, file_format=input_file_format
        )
    except InvalidInputFileFormat:
        return logger.error(f"""Format type could not be read.\nPlease check acceptable input file formats: {", ".join([file_format.value for file_format in FileFormat])}""")
    
    sequence_type = sequence_type or get_seq_type(alignment)

    # set output file format if not specified
    if not output_file_format:
        output_file_format = input_file_format
    else:
        output_file_format = FileFormat[output_file_format]

    # determine smart_gap threshold
    if mode in {TrimmingMode.smart_gap, TrimmingMode.kpi_smart_gap, TrimmingMode.kpic_smart_gap}:
        write_determining_smart_gap_threshold()
        gaps = smart_gap_threshold_determination(alignment, sequence_type, quiet)

    # Print to stdout the user arguments
    write_user_args(
        input_file,
        input_file_format,
        output_file,
        output_file_format,
        sequence_type,
        gaps,
        mode,
        complement,
        use_log,
    )

    # create dictionaries of sequences to keep or trim from the alignment
    keepD, trimD = keep_trim_and_log(
        alignment, gaps, mode, use_log, output_file, complement, sequence_type, quiet
    )

    if use_log:
        warn_if_all_sites_were_trimmed(keepD)

        warn_if_entry_contains_only_gaps(keepD)

    # convert keepD and trimD to multiple sequence alignment objects
    # and write out file
    write_keepD(keepD, output_file, output_file_format)

    # if the -c/--complementary argument was used,
    # create an alignment of the trimmed sequences
    if complement:
        write_trimD(trimD, output_file, output_file_format)

    # print out output statistics
    stats = TrimmingStats(alignment, keepD, trimD)
    write_output_stats(stats, start_time)


def main(argv=None):
    """
    Function that parses and collects arguments              
    """
    # parse and assign arguments
    parser = create_parser()
    args = parser.parse_args()

    # pass to master execute function
    execute(**process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])

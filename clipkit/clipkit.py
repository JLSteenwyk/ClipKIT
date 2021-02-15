#!/usr/bin/env python

import getopt
import logging
import os.path
import sys
import time

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

from .args_processing import process_args
from .files import get_alignment_and_format, FileFormat
from .helpers import (
    keep_trim_and_log,
    write_keepD,
    write_trimD
)
from .modes import TrimmingMode
from .parser import create_parser
from .smart_gap_helper import smart_gap_threshold_determination
from .warnings import (
    check_if_all_sites_were_trimmed,
    check_if_entry_contains_only_gaps,
)
from .write import write_determining_smart_gap_threshold, write_user_args, write_output_stats

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def execute(
    input_file: str,
    input_file_format: FileFormat,
    output_file: str,
    output_file_format: FileFormat,
    gaps: float,
    complement: bool,
    mode: TrimmingMode,
    use_log: bool,
):
    """
    Master execute Function                                      
    This function executes the main functions and calls other    
    subfunctions to trim the input file  
    """
    # logic for whether or not to create a log file
    if use_log:
        # write INFO level logging to file for user
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(f"{output_file}.log", mode="w")
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

    # create start time logger
    start_time = time.time()

    # read in alignment and save the format of the alignment
    alignment, input_file_format = get_alignment_and_format(
        input_file, file_format=input_file_format
    )
    
    # set output file format if not specified
    if not output_file_format:
        output_file_format = input_file_format
    else:
        output_file_format = FileFormat[output_file_format]

    # determine smart_gap threshold
    if mode in {TrimmingMode.smart_gap, TrimmingMode.kpi_smart_gap, TrimmingMode.kpic_smart_gap}:
        write_determining_smart_gap_threshold()
        gaps = smart_gap_threshold_determination(alignment)

    # Print to stdout the user arguments
    write_user_args(
        input_file,
        input_file_format,
        output_file,
        output_file_format,
        gaps,
        mode,
        complement,
        use_log,
    )

    # create dictionaries of sequences to keep or trim from the alignment
    keepD, trimD = keep_trim_and_log(
        alignment, gaps, mode, use_log, output_file, complement
    )

    # check if resulting alingment length is 0
    check_if_all_sites_were_trimmed(keepD)

    # checking if any sequence entry contains only gaps
    check_if_entry_contains_only_gaps(keepD)

    # convert keepD and trimD to multiple sequence alignment objects
    # and write out file
    write_keepD(keepD, output_file, output_file_format)

    # if the -c/--complementary argument was used,
    # create an alignment of the trimmed sequences
    if complement:
        write_trimD(trimD, output_file, output_file_format)

    # print out output statistics
    write_output_stats(alignment, keepD, trimD, start_time)


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

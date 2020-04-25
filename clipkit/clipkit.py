#!/usr/bin/env python

import getopt
import logging
import os.path
import sys
import textwrap
import time

from argparse import ArgumentParser, RawTextHelpFormatter, SUPPRESS, RawDescriptionHelpFormatter
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .helpers import keep_trim_and_log, write_keepD, write_trimD
from .files import get_alignment_and_format, FileFormat
from .modes import TrimmingMode
from .args_processing import process_args

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)


## TODO: Write integration test for log file

####################################################################
### Master execute Function                                      ###
### This function executes the main functions and calls other    ###
### subfunctions to trim the input file  						 ###
####################################################################
def execute(
    inFile: str,
    outFile: str,
    inFileFormat: FileFormat,
    outFileFormat: FileFormat,
    gaps: float,
    complement: bool,
    mode: TrimmingMode,
    use_log: bool,
    ):
    """
    """
    # logic for whether or not to create a log file
    if use_log:
        # write INFO level logging to file for user
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(f"{outFile}.log", mode="w")
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

    # create start time logger
    start_time = time.time()

    # read in alignment and save the format of the alignment
    alignment, inFileFormat = get_alignment_and_format(inFile, file_format=inFileFormat)

    # set output file format if not specified
    if not outFileFormat:
        outFileFormat = inFileFormat

    print(textwrap.dedent(f"""\
        -------------
        | Arguments |
        -------------
        Input file: {inFile} (format: {inFileFormat.value})
        Output file: {outFile} (format: {outFileFormat.value})
        Gaps threshold: {gaps}
        Trimming mode: {mode.value}
        Create complementary output: {complement}
        Create log file: {outFile + '.log' if use_log else False}
    """))

    # create dictionaries of sequences to keep or trim from the alignment
    keepD, trimD = keep_trim_and_log(alignment, gaps, mode, use_log)

    # check if resulting alingment length is 0
    if not len(next(iter(keepD.values()))):
        logger.warning("WARNING: All sites trimmed from alignment. Please use different parameters.")

    # checking if any sequence entry contains only gaps
    for entry, sequence in keepD.items():
        chars_in_sequence = set(sequence)
        if len(chars_in_sequence) == 1 and '-' in chars_in_sequence:
            logger.warning(f"WARNING: header id '{entry}' contains only gaps")
            break
    
    # convert keepD and trimD to multiple sequence alignment objects
    # and write out file
    write_keepD(keepD, outFile, outFileFormat)

    # if the -c/--complementary argument was used,
    # create an alignment of the trimmed sequences
    if complement:
        write_trimD(trimD, outFile, outFileFormat)

    alignment_length = alignment.get_alignment_length()
    output_len  = len(next(iter(keepD.values())))
    trimmed_len = len(next(iter(trimD.values())))
    print(textwrap.dedent(f"""\


        ---------------------
        | Output Statistics |
        ---------------------
        Number of sites kept: {output_len}
        Number of sites trimmed: {trimmed_len}

        Percentage of alignment trimmed: {round((trimmed_len / alignment_length) * 100, 3)}%

        Execution time: {round(time.time() - start_time, 3)}s
    """))



####################################################################
### END Master execute Function 					             ###
####################################################################

####################################################################
### Function that saves arguments and passes them to a           ###
### processing function                                          ###
####################################################################
def main(argv=None):
    """
    Parses arguments 
    """

    parser = ArgumentParser(
        add_help=False,
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        description=textwrap.dedent("""\
              _____ _ _       _  _______ _______ 
             / ____| (_)     | |/ /_   _|__   __|
            | |    | |_ _ __ | ' /  | |    | |   
            | |    | | | '_ \|  <   | |    | |   
            | |____| | | |_) | . \ _| |_   | |   
             \_____|_|_| .__/|_|\_\_____|  |_|   
                       | |                       
                       |_|  

            Citation: Steenwyk et al. Journal, journal info, link

            ClipKIT trims multiple sequence alignments and maintains phylogenetically informative sites.

            Usage: clipkit <input> [optional arguments]
            """)
    )

    # if no arguments are given, print help and exit
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit()

    ## required arguments
    required = parser.add_argument_group("required arguments", description=textwrap.dedent("""\
        <input>                                     input file
                                                    (must be the first argument)
        """)
    )

    required.add_argument(
        "input",
        type=str,
        help=SUPPRESS
    )

    ## optional arguments
    optional = parser.add_argument_group("optional arguments", description=textwrap.dedent("""\
        -o, --output <output_file_name>             output file name 
                                                    (default: input file named with '.clipkit' suffix)

        -m, --modes <kpi, gappy, kpi-gappy>         trimming mode 
                                                    (default: gappy)

        -g, --gaps <threshold of gaps>              specifies gaps threshold
                                                    (default: 0.9)

        -if, --input_file_format <file_format>      specifies input file format
                                                    (default: auto-detect)    

        -of, --output_file_format <file_format>     specifies output file format
                                                    (default: same as input file format)

        -l, --log                                   creates a log file
                                                    (input file named with '.log' suffix)

        -c, --complementary                         creates complementary alignment of trimmed sequences
                                                    (input file named with '.log' suffix)
 
        -h, --help                                  help message
        -v, --version                               print version

        -------------------------------------
        | Detailed explanation of arguments | 
        -------------------------------------
        Modes
            kpi: keep only parsimony informative sites in an alignment
            gappy: trim sites that are greater than the gaps threshold
            kpi-gappy: a combination of kpi- and gappy-based trimming

        Gaps
            Positions with gappyness greater than threshold will be trimmed. 
            Must be between 0 and 1. (Default: 0.9). This argument is ignored
            when using the kpi mode of trimming.

        Input and output file formats
            Supported input and output files include:
            fasta, clustal, maf, mauve, phylip, phylip-sequential, 
            phylip-relaxed, and stockholm

        Log
            Creates a log file that summarizes the characteristics of each position.
            The log file has four columns.
            - Column 1 is the position in the alignment (starting at 1), 
            - Column 2 reports if the site was trimmed or kept (trim and keep, respectively),
            - Column 3 reports if the site is a parsimony informative site or not (PI and nPI, respectively), and
            - Column 4 reports the gappyness of the the position (number of gaps / entries in alignment)

        Complementary
            Creates an alignment file of only the trimmed sequences
        """)
    )

    optional.add_argument(
        "-o",
        "--output",
        help=SUPPRESS,
        metavar="output"
    )

    optional.add_argument(
        "-m",
        "--mode",
        help=SUPPRESS,
        nargs="?",
        choices=("kpi", "gappy", "kpi-gappy")
    )

    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help=SUPPRESS,
    )

    optional.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 0.1",
        help=SUPPRESS,
    )

    optional.add_argument(
        "-g",
        "--gaps",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="fraction of gaps",
    )

    file_format_choices = [file_format.value.lower() for file_format in FileFormat]
    optional.add_argument(
        "-if",
        "--input_file_format",
        type=str,
        required=False,
        choices=file_format_choices,
        help=SUPPRESS,
        metavar=""
    )

    optional.add_argument(
        "-of",
        "--output_file_format",
        type=str,
        required=False,
        choices=file_format_choices,
        help=SUPPRESS,
        metavar=''
    )

    optional.add_argument(
        "-l",
        "--log",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "-c",
        "--complementary",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    # parse and assign arguments
    args = parser.parse_args()

    # pass to master execute function
    execute(*process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])

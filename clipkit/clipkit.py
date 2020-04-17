#!/usr/bin/env python

import logging

import sys
import getopt
import os.path
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from argparse import ArgumentParser, RawTextHelpFormatter
from .helpers import keep_trim_and_log, write_keepD, write_trimD
from .files import automatic_file_type_determination, FileFormat

from .modes import TrimmingMode


ck_log = logging.getLogger(__name__)

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
    use_log: bool,
    gaps: float,
    complement: bool,
    mode: TrimmingMode
    ):
    """
    """

    # read in alignment and save the format of the alignment

    ## TODO: Jacob -- Refactor with new argument to automatic_file_type_determination
    ## def automatic_file_type_determination(inFile, file_format=None)
    ## call it here with automatic_file_type_determination(inFile, file_format=inFileFormat)
    ## rename function to also note that the alignment is read in -- e.g., get_alignment_and_format()
    if inFileFormat:
        alignment = AlignIO.read(open(inFile), inFileFormat)
    else:
        alignment, inFileFormat = automatic_file_type_determination(inFile)
    
    # set output file format if not specified
    if not outFileFormat:
        outFileFormat = inFileFormat

    # create dictionaries of sequences to keep or trim from the alignment
    keepD, trimD, logArr = keep_trim_and_log(alignment, gaps, mode)

    # convert keepD and trimD to multiple sequence alignment objects
    # and write out file
    write_keepD(keepD, outFile, outFileFormat)
    
    # if the -c/--complementary argument was used, 
    # create an alignment of the trimmed sequences
    if complement:
        write_trimD(trimD, outFile, outFileFormat)

    # if the -l/--log argument was used,
    # create a log file with information about each
    # position in the alignment
    # TODO: handle log output with Python logger
    if use_log:
        outFileLog = outFile + ".log"
        with open(outFileLog, 'w+') as f:
            f.writelines('\t'.join(entry)+'\n' for entry in logArr)

####################################################################
### END Master execute Function 					             ###
####################################################################

####################################################################
### Function that read the input files and create output         ###
### Input files include gene trees and species tree              ###
### Output is printed to stdout                                  ###
####################################################################
def main(
    argv=None
    ):
    """
    Parses arguments 
    """
    if not argv:
        argv = sys.argv[1:]

    parser = ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter,
    	usage="""python %(prog)s -i input -o output
    
    optional arguments
    ------------------
    [-g gap_threshold] [-if input file format] [-m kpi, gappy, or kpi-gappy]
    [-of output file format] [-l log file] [-c complementary sequence]""",
    description="""Citation: Steenwyk et al. Journal, journal info, link

    Fast trimming of multiple sequence alignments to maintain phylogenetically informative sites. ClipKIT 
    will...""")

    ## required arguments
    # create a required group of arguments
    required = parser.add_argument_group('required arguments')
    # required input and output arguments
    required.add_argument('-i', '--input', required=True, 
        help="""Input file. See file_format for acceptable file formats.
        """)

    ## optional arguments
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument('-o', '--output', help="Output file. Output format will match input format.")
    
    optional.add_argument('-m', '--mode', help='kpi, gappy, or kpi-gappy', nargs='?', 
        choices=('kpi', 'gappy', 'kpi-gappy'))
    
    optional.add_argument('-h', '--help', action='help',
        help="""show this help message and exit
        """)

    optional.add_argument('-v', '--version', action='version', 
        version='%(prog)s 0.1', help="""print program version
    	""")
    
    optional.add_argument('-g', '--gaps', type=float, required=False,
        help="""Specifies gaps threshold. Must be between 0 and 1. (Default: 0.9)
        """)

    file_format_choices = [file_format.value.lower() for file_format in FileFormat]
    optional.add_argument('-if', '--input_file_format', type=str, required=False,
        choices=file_format_choices,
        help="""Specify input file format. Currently accepted file formats are: 
    fasta, clustal, maf, mauve, phylip, phylip-sequential, phylip-relaxed, stockholm
    (Default: auto-detect) 
    """)

    optional.add_argument('-of', '--output_file_format', type=str, required=False,
        choices=file_format_choices,
        help="""Specify output file format. Currently accepted file formats are: 
    fasta, clustal, maf, mauve, phylip, phylip-sequential, phylip-relaxed, stockholm
    (Default: input file type)
    """)

    optional.add_argument('-l', '--log', action='store_true', required=False,
        help="""Creates a log file that summarizes each position (suffix = '.log').
    The log file has four columns. Column 1 is the position in the alignment (starting at 1), 
    column 2 reports if the site was trimmed or kept (trim and keep, respectively), column 3
    reports if the site is a parsimony informative site or not (PI and nPI, respectively),
    and column 4 reports the gappyness of the the position (number of gaps / entries in alignment)
    """)

    optional.add_argument('-c', '--complementary', action='store_true', required=False, 
        help="""Creates an alignment file of the trimmed sequences (suffix = '.complement')
        """)
    
    # parse and assign arguments
    args = parser.parse_args()
    inFile = args.input
    outFile = args.output or f"{inFile}.clipkit"

    if inFile == outFile:
        ck_log.info("Input and output files can't have the same name.")
        sys.exit()

    ## assign optional arguments
    mode = args.mode or TrimmingMode.gappy
    use_log = args.log or False
    gaps = float(args.gaps) if args.gaps else 0.9
    complement = args.complementary or False
    
    inFileFormat = args.input_file_format
    outFileFormat = args.output_file_format

    # pass to master execute function
    execute(
        inFile,
        outFile,
        inFileFormat,
        outFileFormat,
        use_log,
        gaps,
        complement,
        mode
        )

if __name__ == '__main__':
    main(sys.argv[1:])

# clipkit input_path
# clipkit input_path output_path
# cliplit input_path -o output_path
# cliplit -i input_path -o output_path
# cliplit ... -m mode , ..., 

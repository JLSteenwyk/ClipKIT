#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from argparse import ArgumentParser, RawTextHelpFormatter

from .helpers import automatic_file_type_determination, keep_trim_and_log

####################################################################
### Master execute Function                                      ###
### This function executes the main functions and calls other    ###
### subfunctions to trim the input file  						 ###
####################################################################
def execute(
    inFile,
    outFile,
    log,
    gaps,
    complement,
    inFileFormat,
    outFileFormat,
    mode
    ):
    """
    """

    # read in alignment and save the format of the alignment
    if inFileFormat != 'NA':
        alignment = AlignIO.read(open(inFile), inFileFormat)
    elif inFileFormat == 'NA':
        alignment, inFileFormat = automatic_file_type_determination(inFile)
    
    # set output file format if not specified
    if outFileFormat == 'NA':
        outFileFormat = inFileFormat

    # create dictionaries of sequences to keep or trim from the alignment
    keepD, trimD, logArr = keep_trim_and_log(alignment, gaps, mode)

    # convert keepD and trimD to multiple sequence alignment objects
    # and write out file
    seqList = []
    for indiv, seq in keepD.items():
        seqList.append(SeqRecord(Seq(str(keepD[indiv])), id=str(indiv), description=''))
    keepMSA = MultipleSeqAlignment(seqList)
    SeqIO.write(keepMSA, outFile, outFileFormat)
    
    # if the -c/--complementary argument was used, 
    # create an alignment of the trimmed sequences
    if complement:
        seqList = []
        for indiv, seq in keepD.items():
            seqList.append(SeqRecord(Seq(str(trimD[indiv])), id=str(indiv), description=''))
        trimMSA = MultipleSeqAlignment(seqList)
        completmentOut = str(outFile) + ".complement"
        SeqIO.write(trimMSA, completmentOut, outFileFormat)

    # if the -l/--log argument was used,
    # create a log file with information about each
    # position in the alignment
    if log:
        outFileLog = outFile + ".log"
        with open(outFileLog, 'w+') as f:
            f.writelines('\t'.join(entry)+'\n' for entry in logArr)

####################################################################
### END Master execute Function 					             ###
####################################################################


## Print help message if an invalid string was used to specify
## the input or output file format
def help_wrong_file_format(
	fileFormat
    ):
    """
    Prints a message that the wrong file format was specified
    and then prints out a list of acceptable file formats

    Arguments
    ---------
    argv: fileFormat
        user specified input or output file format  
    """

    print("\n")
    print(fileFormat, "is not an accepted alignment format. Accepted alignment formats are:")
    print("'fasta', 'clustal', 'maf', 'mauve', 'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm'\n")
    print("For more help, use the -h/--help argument\n")
    

####################################################################
### Functions that read the input files and create output        ###
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
        # TODO: clean up
        argv = sys.argv[1:]

    # initialize argument variables
    inFile        = ''
    outFile       = ''
    log           = ''
    gaps          = ''
    complement    = ''
    inFileFormat  = ''
    outFileFormat = ''

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

    required.add_argument('-o', '--output', required=True,
        help="""Output file. Output format will match input format.
        """)

    ## optional arguments
    optional = parser.add_argument_group('optional arguments')
    
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

    optional.add_argument('-if', '--input_file_format', required=False,
        help="""Specify input file format. Currently accepted file formats are: 
    fasta, clustal, maf, mauve, phylip, phylip-sequential, phylip-relaxed, stockholm
    (Default: auto-detect) 
    """)

    optional.add_argument('-of', '--output_file_format', required=False,
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
    outFile = args.output

    if inFile == outFile:
        print("Input and output files can't have the same name.")
        sys.exit()

    ## assign optional arguments
    if args.mode:
        mode = args.mode
    else:
        mode = 'gappy'
    if args.log:
        log = args.log
    else:
        log = False
    
    if args.gaps:
        gaps = float(args.gaps)
    else:
        gaps = 0.9
    
    if args.complementary:
        complement = args.complementary
    else:
        complement = False
    
    fileFormats = ['fasta', 'clustal', 'maf', 'mauve', 'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm']
    if args.input_file_format:
        inFileFormat = args.input_file_format
        
        if inFileFormat not in fileFormats:
            help_wrong_file_format(inFileFormat)
            sys.exit()

    else:
        inFileFormat = 'NA'

    if args.output_file_format:
        outFileFormat = args.output_file_format
        
        if outFileFormat not in fileFormats:
            help_wrong_file_format(outFileFormat)
            sys.exit()
    else:
        outFileFormat = 'NA'

    # pass to master execute function
    execute(
        inFile,
        outFile,
        log,
        gaps,
        complement,
        inFileFormat,
        outFileFormat,
        mode
        )

if __name__ == '__main__':
    main(sys.argv[1:])
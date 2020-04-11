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


####################################################################
### Supporting Functions 	                                     ###
### This block of code contains all supporting functions that    ###
### read the input files, calculate gappyness, etc  			 ###
####################################################################
## Function to automatically determine the format of the alignment file
## and read in the alignment. Returns alignment object and fileFormat
def automatic_file_type_determination(
    inFile
    ):
    """
    Automatically determines what type of input file was used
    and reads in the alignment file

    Arguments
    ---------
    argv: inFile
        input file specified with -i, --input
    """

    # save list of different file formats
    fileFormats = ['fasta', 'clustal', 'maf', 'mauve',
        'phylip', 'phylip-sequential', 'phylip-relaxed', 'stockholm']

    # loop through file formats and attempt to read in file in that format
    for fileFormat in fileFormats:
        try:
            alignment = AlignIO.read(open(inFile), fileFormat)
            return(alignment, fileFormat)
            break
        # the following exceptions refer to skipping over errors 
        # associated with reading the wrong input file
        except ValueError:
            continue
        except AssertionError:
            continue
   
## Function to determine which positions of an alignment should be 
## kept or trimmed 
def keep_trim_and_log(
    alignment, gaps, mode
    ):
    """
    Determines positions to keep or trim and saves these positions
    to dictionaries named keepD and trimD. For both dictionaries,
    keys are the names of the sequence entries and the values are
    sequences that will be kept or trimmed

    Arguments
    ---------
    argv: alignment
        biopython multiple sequence alignment object
    argv: gaps
        gaps threshold to determine if a position is too gappy or not
    """

    # initialize dictionaries that will eventually be populated with
    # alignment positions to keep or trim (keys) and the sequence at
    # that position (values). Also, initialize a list of log information
    # that will be kept in an array format
    keepD = {}
    for entry in alignment:
        keepD[entry.id] = []
    trimD = {}
    for entry in alignment:
        trimD[entry.id] = []
    logArr = []

    # loop through alignment
    for i in range(0, alignment.get_alignment_length(), int(1)):
        
        # save the sequence at the position to a string 
        seqAtPosition  = ''
        seqAtPosition += alignment[:, i]
        seqAtPosition  = seqAtPosition.upper()

        # determine the length and number of gaps in an alignment position
        lengthOfSeq = len(seqAtPosition)
        numOfGaps   = seqAtPosition.count('-')
        gappyness   = numOfGaps/lengthOfSeq


        ## determine if the site is parsimony informative 
        ## "A site is parsimony-informative if it contains at least two types of nucleotides 
        ## (or amino acids), and at least two of them occur with a minimum frequency of two."
        ## https://www.megasoftware.net/web_help_7/rh_parsimony_informative_site.htm

        # Create a dictionary that tracks the number of occurences of each character 
        # excluding gaps or '-'
        numOccurences = {}
        for char in set(seqAtPosition.replace('-','')):	
            numOccurences[char] = seqAtPosition.count(char)

        # if the number of values that are greater than two 
        # in the numOccurences dictionary is greater than two, 
        # the site is parsimony informative
        d = dict((k, v) for k, v in numOccurences.items() if v >= 2)
        if len(d) >= 2:
            parsimony_informative = True
        else:
            parsimony_informative = False

        if mode == 'kpi-gappy':
            # if gappyness is lower than gaps threshold and the site is
            # parismony informative, save to keepD. Otherwise, save to trimD.
            # All the while, save information to logD
            if gappyness <= gaps and parsimony_informative:
                # save to keepD
                for entry in alignment:
                    keepD[entry.id].append(entry.seq._data[i])
                # save to logL - structure is site, parsimony informative (PI) or not (nPI)
                # gappyness, kept or trimmed
                temp = []
                temp.append(str(i+1))
                temp.append('keep')
                temp.append("PI")
                temp.append(str(gappyness))
                logArr.append(temp)
            else:
                # save to trimD 
                for entry in alignment:
                    trimD[entry.id].append(entry.seq._data[i])
                # save to logL - structure is site, parsimony informative (PI) or not (nPI)
                # gappyness, kept or trimmed
                temp = []
                temp.append(str(i+1))
                temp.append('trim')
                temp.append("nPI")
                temp.append(str(gappyness))
                logArr.append(temp)
        elif mode == 'gappy':
            # if gappyness is lower than gaps threshold and the site is
            # parismony informative, save to keepD. Otherwise, save to trimD.
            # All the while, save information to logD
            if gappyness <= gaps:
                # save to keepD
                for entry in alignment:
                    keepD[entry.id].append(entry.seq._data[i])
                # save to logL - structure is site, parsimony informative (PI) or not (nPI)
                # gappyness, kept or trimmed
                if parsimony_informative:
                    temp = []
                    temp.append(str(i+1))
                    temp.append('keep')
                    temp.append("PI")
                    temp.append(str(gappyness))
                    logArr.append(temp)
                else: 
                    temp = []
                    temp.append(str(i+1))
                    temp.append('keep')
                    temp.append("nPI")
                    temp.append(str(gappyness))
                    logArr.append(temp)
            else:
                # save to trimD 
                for entry in alignment:
                    trimD[entry.id].append(entry.seq._data[i])
                # save to logL - structure is site, parsimony informative (PI) or not (nPI)
                # gappyness, kept or trimmed
                if parsimony_informative:
                    temp = []
                    temp.append(str(i+1))
                    temp.append('trim')
                    temp.append("PI")
                    temp.append(str(gappyness))
                    logArr.append(temp)
                else: 
                    temp = []
                    temp.append(str(i+1))
                    temp.append('trim')
                    temp.append("nPI")
                    temp.append(str(gappyness))
                    logArr.append(temp)
        
        elif mode == 'kpi':
            # if site is parismony informative, save to keepD. 
            # Otherwise, save to trimD.
            # All the while, save information to logD
            if parsimony_informative:
                # save to keepD
                for entry in alignment:
                    keepD[entry.id].append(entry.seq._data[i])
 
                temp = []
                temp.append(str(i+1))
                temp.append('keep')
                temp.append("PI")
                temp.append(str(gappyness))
                logArr.append(temp)

            else:
                # save to trimD 
                for entry in alignment:
                    trimD[entry.id].append(entry.seq._data[i])

                temp = []
                temp.append(str(i+1))
                temp.append('trim')
                temp.append("nPI")
                temp.append(str(gappyness))
                logArr.append(temp)
            
    # join elements in value lists in keepD and trimD 
    for k, v in keepD.items():
        keepD[k] = ''.join(v)
    for k, v in trimD.items():
        trimD[k] = ''.join(v)

    return(keepD, trimD, logArr)

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
    argv
    ):
    """
    Parses arguments 
    """

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

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from .modes import kpi_gappy_mode
from .modes import gappy_mode
from .modes import kpi_mode
from .files import FileFormat

####################################################################
### Supporting Functions 	                                     ###
### This block of code contains all supporting functions that    ###
### read the input files, calculate gappyness, etc  			 ###
####################################################################

## Function to get the sequence at a given column. Function
## will also determine features of the position including 
## length and number of gaps.
## Support function for keep_trim_and_log()
def get_sequence_at_position_and_report_features(
    alignment, i
    ):
    """
    Count the occurence of each character at a given position
    in an alignment. This information is used to determine
    if the aligment is parsimony informative or not. When 
    count characters, gaps ('-') are excluded

    Arguments
    ---------
    argv: seqAtPosition
        string that contains the sequence at a given column 
    """

    # save the sequence at the position to a string 
    seqAtPosition  = ''
    seqAtPosition += alignment[:, i]
    seqAtPosition  = seqAtPosition.upper()

    # determine the length and number of gaps in an alignment position
    lengthOfSeq = len(seqAtPosition)
    numOfGaps   = seqAtPosition.count('-')
    gappyness   = numOfGaps/lengthOfSeq

    return seqAtPosition, gappyness


## Function to count the number of occurences in each character
## Support function for keep_trim_and_log()
def count_characters_at_position(
    seqAtPosition
    ):
    """
    Count the occurence of each character at a given position
    in an alignment. This information is used to determine
    if the aligment is parsimony informative or not. When 
    count characters, gaps ('-') are excluded

    Arguments
    ---------
    argv: seqAtPosition
        string that contains the sequence at a given column
    """

    numOccurences = {}
    for char in set(seqAtPosition.replace('-','')): 
        numOccurences[char] = seqAtPosition.count(char)
    return numOccurences


## Function determine if a site is parsimony informative or not
## Support function for keep_trim_and_log()
def determine_if_parsimony_informative(
    numOccurences
    ):
    """
    Determines if a site is parsimony informative.
    A site is parsimony-informative if it contains at least two types of nucleotides 
    (or amino acids), and at least two of them occur with a minimum frequency of two."
    https://www.megasoftware.net/web_help_7/rh_parsimony_informative_site.htm

    Arguments
    ---------
    argv: numOccurences
        dictionary with sequence characters (keys) and their counts (values)
    """

    d = dict((k, v) for k, v in numOccurences.items() if v >= 2)
    if len(d) >= 2:
        parsimony_informative = True
    else:
        parsimony_informative = False

    return parsimony_informative

## Function to initialize and populate keepD and trimD with entry ids  
## and empty list elements. This will eventually contain the 
## sequences that are trimmed or kept by clipkit. Lastly, initialize
## an array to keep log information
def populate_empty_keepD_and_trimD(
    alignment
    ):
    """
    Creates barebones dictionaries for sites kept and trimmed. Creates
    an array for keeping log information.

    Arguments
    ---------
    argv: alignment
        biopython multiple sequence alignment object
    """
    keepD = {}
    for entry in alignment:
        keepD[entry.id] = []
    trimD = {}
    for entry in alignment:
        trimD[entry.id] = []
    logArr = []

    return keepD, trimD, logArr

## Function to join the elements of keepD and trimD into a nicer
## arrangement of key value pairs
def join_keepD_and_trimD(
    keepD, trimD
    ):
    """
    Currently, each position is its own element. This function
    will join those elements into one string.

    Arguments
    ---------
    argv: keepD
        dictionary of sequences to be kept after trimmed
    argv: trimD
        dictionary of sequences to be trimmed off
    """

    # join elements in value lists in keepD and trimD 
    for k, v in keepD.items():
        keepD[k] = ''.join(v)
    for k, v in trimD.items():
        trimD[k] = ''.join(v)

    return keepD, trimD

## Function to write out keepD to output file
# TODO: write unit test
def write_keepD(
    keepD,
    outFile,
    outFileFormat: FileFormat
    ):
    """
    This creates a biopython multisequence alignment object. Object
    is populated with sites that are kept after trimming is finished

    Arguments
    ---------
    argv: keepD
        dictionary of sites that are kept in resulting alignment
    argv: outFile
        name of the output file
    argv: outFileFormat
        output file format
    """

    seqList = []
    for indiv, seq in keepD.items():
        seqList.append(SeqRecord(Seq(str(keepD[indiv])), id=str(indiv), description=''))
    keepMSA = MultipleSeqAlignment(seqList)
    SeqIO.write(keepMSA, outFile, outFileFormat.value)

## Function to write out trimD to output file
# TODO: write unit test
def write_trimD(
    trimD, completmentOut, outFileFormat
    ):
    """
    This creates a biopython multisequence alignment object. Object
    is populated with sites that are trimmed after trimming is finished

    Arguments
    ---------
    argv: trimD
        dictionary of sites that are trimmed in final alignment
    argv: complementOut
        name of complementary output file
    argv: outFileFormat
        file format of complementary output file
    """

    seqList = []
    for indiv, seq in trimD.items():
        seqList.append(SeqRecord(Seq(str(trimD[indiv])), id=str(indiv), description=''))
    trimMSA = MultipleSeqAlignment(seqList)
    completmentOut = str(outFile) + ".complement"
    SeqIO.write(trimMSA, completmentOut, outFileFormat)

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
    argv: mode
        mode of how to trim alignment
    """

    # initialize dictionaries that will eventually be populated with
    # alignment positions to keep or trim (keys) and the sequence at
    # that position (values). Also, initialize a list of log information
    # that will be kept in an array format
    keepD, trimD, logArr = populate_empty_keepD_and_trimD(alignment)

    # loop through alignment
    for i in range(0, alignment.get_alignment_length(), int(1)):
        
        # save the sequence at the position to a string and calculate the gappyness of the site
        seqAtPosition, gappyness = get_sequence_at_position_and_report_features(alignment, i)

        ## determine if the site is parsimony informative and trim accordingly
        # Create a dictionary that tracks the number of occurences of each character 
        # excluding gaps or '-'
        numOccurences = count_characters_at_position(seqAtPosition)

        # if the number of values that are greater than two 
        # in the numOccurences dictionary is greater than two, 
        # the site is parsimony informative
        parsimony_informative = determine_if_parsimony_informative(numOccurences)

        # depending on the mode, trim the alignment
        if mode == 'kpi-gappy':
            keepD, trimD, logArr = kpi_gappy_mode(
                gappyness, parsimony_informative, 
                keepD, trimD, logArr, i, gaps, alignment
                )
        elif mode == 'gappy':
            keepD, trimD, logArr = gappy_mode(
                gappyness, parsimony_informative, 
                keepD, trimD, logArr, i, gaps, alignment
                )
        elif mode == 'kpi':
            keepD, trimD, logArr = kpi_mode(
                gappyness, parsimony_informative,
                keepD, trimD, logArr, i, gaps, alignment
                )

    # join elements in value lists in keepD and trimD 
    keepD, trimD = join_keepD_and_trimD(keepD, trimD)

    return keepD, trimD, logArr
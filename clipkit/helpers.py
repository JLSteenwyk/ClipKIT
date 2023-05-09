import re
import sys
import textwrap
import time

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from math import floor
from tqdm import tqdm

from .modes import TrimmingMode, trim
from .files import FileFormat
from .write import write_processing_aln, write_output_files_message

from enum import Enum


class SeqType(Enum):
    aa = "aa"
    nt = "nt"

class SiteClassificationType(Enum):
    parsimony_informative = "parsimony-informative"
    constant = "constant"
    singleton = "singleton"
    other = "other"

def remove_gaps(seq: str) -> str:
    return re.sub("\-|\?|\*|x", "", seq, flags=re.IGNORECASE)

def get_seq_type(
    alignment: MultipleSeqAlignment,
) -> SeqType:
    seq = str(alignment[0].seq)
    seq = remove_gaps(seq)
    if len(seq) < 200:
        seq = "".join([str(record.seq) for record in alignment])
        seq = remove_gaps(seq)

    if len(set(seq.upper())) > 5:
        sequence_type = SeqType.aa
    else:
        sequence_type = SeqType.nt

    return sequence_type


def get_sequence_at_position_and_report_features(alignment, i, char):
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
    seqAtPosition = ""
    seqAtPosition += alignment[:, i]
    seqAtPosition = seqAtPosition.upper()

    # determine the length and number of gaps in an alignment position
    lengthOfSeq = len(seqAtPosition)
    if char == SeqType.aa:
        for gap_char in ["?", "*", "X"]:
            seqAtPosition = seqAtPosition.replace(gap_char, "-")
    else:
        for gap_char in ["?", "*", "X", "N"]:
            seqAtPosition = seqAtPosition.replace(gap_char, "-")

    numOfGaps = seqAtPosition.count("-")
    gappyness = numOfGaps / lengthOfSeq

    return seqAtPosition, gappyness


def count_characters_at_position(seqAtPosition):
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
    for char in set(seqAtPosition.replace("-", "")):
        numOccurences[char] = seqAtPosition.count(char)
    return numOccurences


def determine_site_classification_type(numOccurences) -> SiteClassificationType:
    """
    Determines if a site is parsimony informative or constant.
    A site is parsimony-informative if it contains at least two types of nucleotides 
    (or amino acids), and at least two of them occur with a minimum frequency of two.
    https://www.megasoftware.net/web_help_7/rh_parsimony_informative_site.htm

    A site is constant if it contains only one character and that character occurs
    at least twice. https://www.megasoftware.net/web_help_7/rh_constant_site.htm

    Arguments
    ---------
    argv: numOccurences
        dictionary with sequence characters (keys) and their counts (values)
    """
    # create a dictionary of characters that occur at least twice
    d = dict((k, v) for k, v in numOccurences.items() if v >= 2)
    if len(d) >= 2:
        return SiteClassificationType.parsimony_informative
    elif len(d) == 1 and len(numOccurences) == 1:
        return SiteClassificationType.constant
    elif len(d) == 1 and len(numOccurences) > 1:
        return SiteClassificationType.singleton

    return SiteClassificationType.other


def populate_empty_keepD_and_trimD(alignment):
    """
    Creates barebones dictionaries for sites kept and trimmed. Creates
    an array for keeping log information.

    Arguments
    ---------
    argv: alignment
        biopython multiple sequence alignment object
    """
    keepD = {}
    trimD = {}
    alignment_length = alignment.get_alignment_length()
    for entry in alignment:
        keepD[entry.description] = np.zeros([alignment_length], dtype=bytes)
        trimD[entry.description] = np.zeros([alignment_length], dtype=bytes)

    return keepD, trimD


def join_keepD_and_trimD(keepD, trimD):
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
        keepD[k] = "".join(np.char.decode(v))
    for k, v in trimD.items():
        trimD[k] = "".join(np.char.decode(v))

    return keepD, trimD


def write_keepD(keepD, outFile, outFileFormat: FileFormat):
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
    for indiv in keepD.keys():
        seqList.append(SeqRecord(Seq(str(keepD[indiv])), id=str(indiv), description=""))
    keepMSA = MultipleSeqAlignment(seqList)
    if outFileFormat.value == 'phylip_relaxed':
        SeqIO.write(keepMSA, outFile, 'phylip-relaxed')
    elif outFileFormat.value == 'phylip_sequential':
        SeqIO.write(keepMSA, outFile, 'phylip-sequential')
    else:
        SeqIO.write(keepMSA, outFile, outFileFormat.value)


def write_trimD(trimD, outFile: str, outFileFormat: FileFormat):
    """
    This creates a biopython multisequence alignment object. Object
    is populated with sites that are trimmed after trimming is finished

    Arguments
    ---------
    argv: trimD
        dictionary of sites that are trimmed in final alignment
    argv: outFileFormat
        file format of complementary output file
    argv: name of output file
    """

    seqList = []
    for indiv in trimD.keys():
        seqList.append(SeqRecord(Seq(str(trimD[indiv])), id=str(indiv), description=""))
    trimMSA = MultipleSeqAlignment(seqList)
    completmentOut = str(outFile) + ".complement"
    SeqIO.write(trimMSA, completmentOut, outFileFormat.value)


def keep_trim_and_log(
    alignment, gaps: float, mode: TrimmingMode, use_log: bool, outFile, complement, char: SeqType, quiet: bool
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
    keepD, trimD = populate_empty_keepD_and_trimD(alignment)

    alignment_length = alignment.get_alignment_length()

    # loop through alignment
    write_processing_aln()
    for i in tqdm(range(alignment_length), disable=quiet, postfix="trimmer"):
        # save the sequence at the position to a string and calculate the gappyness of the site
        seqAtPosition, gappyness = get_sequence_at_position_and_report_features(
            alignment, i, char
        )

        ## determine if the site is parsimony informative and trim accordingly
        # Create a dictionary that tracks the number of occurences of each character
        # excluding gaps or '-'
        numOccurences = count_characters_at_position(seqAtPosition)

        # if the number of values that are greater than two
        # in the numOccurences dictionary is greater than two,
        # the site is parsimony informative
        #
        # determine if a site is a constant site or not. A constant site is
        # defined as a site that contains the same nucl or amino acid at all
        # sequence entries. Additionally, that nucl or amino acid must occur at
        # least twice
        site_classification_type = determine_site_classification_type(
            numOccurences
        )

        # trim based on the mode
        keepD, trimD = trim(
            gappyness,
            site_classification_type,
            keepD,
            trimD,
            i,
            gaps,
            alignment,
            mode,
            use_log,
        )

    # print to stdout that output files are being written
    write_output_files_message(outFile, complement, use_log)
    # join elements in value lists in keepD and trimD
    keepD, trimD = join_keepD_and_trimD(keepD, trimD)

    return keepD, trimD

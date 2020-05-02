import logging
from enum import Enum
from typing import Tuple

from Bio import AlignIO

log = logging.getLogger(__name__)


class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    phylip_rel = "phylip-relaxed"
    stockholm = "stockholm"


def get_alignment_and_format(inFile: str, file_format: FileFormat):
    """
    Automatically determines what type of input file was used
    and reads in the alignment file

    Arguments
    ---------
    argv: inFile
        input file specified with -i, --input
    """

    # if file format is provided, read the file according to the user's file format
    if file_format:
        alignment = AlignIO.read(open(inFile), file_format.value)
        return alignment, file_format
    else:
        # loop through file formats and attempt to read in file in that format
        for fileFormat in FileFormat:
            try:
                alignment = AlignIO.read(open(inFile), fileFormat.value)
                return alignment, fileFormat
            # the following exceptions refer to skipping over errors
            # associated with reading the wrong input file
            except ValueError:
                continue
            except AssertionError:
                continue

        raise Exception("Input file could not be read")

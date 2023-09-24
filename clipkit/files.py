from enum import Enum
from .logger import log_file_logger

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from .exceptions import InvalidInputFileFormat


class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_sequential = "phylip_sequential"
    phylip_relaxed = "phylip_relaxed"
    stockholm = "stockholm"


def get_alignment_and_format(
    input_file_name: str, file_format: FileFormat
) -> tuple[MultipleSeqAlignment, FileFormat]:
    """
    Automatically determines what type of input file was used
    and reads in the alignment file
    """

    if file_format:
        file_format = FileFormat(file_format)
        alignment = AlignIO.read(open(input_file_name), file_format.value)
        return alignment, file_format
    else:
        # attempt to auto-detect file format
        for fileFormat in FileFormat:
            try:
                alignment = AlignIO.read(open(input_file_name), fileFormat.value)
                return alignment, fileFormat
            # the following exceptions refer to skipping over errors
            # associated with reading the wrong input file
            except ValueError:
                continue
            except AssertionError:
                continue

        raise InvalidInputFileFormat("File could not be read")


def write_debug_log_file(msa):
    for info in msa.generate_debug_log_info():
        log_file_logger.debug(f"{str(info[0] + 1)} {info[1]} {info[2].value} {info[3]}")

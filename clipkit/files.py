from enum import Enum
from .logger import log_file_logger

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np

from .exceptions import InvalidInputFileFormat
from .ecomp import read_ecomp
from .ecomp.reader import HEADER_MAGIC, EcompDecodeError


class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_sequential = "phylip_sequential"
    phylip_relaxed = "phylip_relaxed"
    stockholm = "stockholm"
    ecomp = "ecomp"


def get_alignment_and_format(
    input_file_name: str, file_format: FileFormat
) -> tuple[MultipleSeqAlignment, FileFormat]:
    """
    Automatically determines what type of input file was used
    and reads in the alignment file
    """

    if file_format:
        file_format = FileFormat(file_format)
        if file_format == FileFormat.ecomp:
            alignment, _ = read_ecomp(input_file_name)
            return alignment, file_format
        with open(input_file_name) as handle:
            alignment = AlignIO.read(handle, file_format.value)
        return alignment, file_format
    else:
        # attempt to auto-detect file format
        if _looks_like_ecomp(input_file_name):
            try:
                alignment, _ = read_ecomp(input_file_name)
                return alignment, FileFormat.ecomp
            except (EcompDecodeError, OSError):
                pass
        for fileFormat in FileFormat:
            if fileFormat == FileFormat.ecomp:
                continue
            try:
                with open(input_file_name) as handle:
                    alignment = AlignIO.read(handle, fileFormat.value)
                return alignment, fileFormat
            # the following exceptions refer to skipping over errors
            # associated with reading the wrong input file
            except ValueError:
                continue
            except AssertionError:
                continue

        raise InvalidInputFileFormat("File could not be read")


def get_custom_sites_to_trim(file_path: str, aln_length: int) -> list:
    with open(file_path) as f:
        lines = f.read().splitlines()

    sites_to_trim = []
    sites_to_keep = []
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            continue

        site = line.split("\t")
        if len(site) != 2:
            raise ValueError(
                f"Invalid CST format on line {line_number}: expected '<position>\\t<trim|keep>'"
            )

        try:
            pos = int(site[0]) - 1
        except ValueError as exc:
            raise ValueError(
                f"Invalid CST position on line {line_number}: '{site[0]}'"
            ) from exc

        if pos < 0 or pos >= aln_length:
            raise ValueError(
                f"CST position out of range on line {line_number}: {pos + 1} (alignment length: {aln_length})"
            )

        action = site[1].strip().lower()
        if action == "trim":
            sites_to_trim.append(pos)
        elif action == "keep":
            sites_to_keep.append(pos)
        else:
            raise ValueError(
                f"Invalid CST action on line {line_number}: '{site[1]}'. Expected 'trim' or 'keep'."
            )

    if len(sites_to_trim) == 0:
        # we only had keeps so treat every other site as a trim
        sites_to_trim = list(
            np.setdiff1d(np.arange(aln_length), np.array(sites_to_keep))
        )

    return sites_to_trim


def write_debug_log_file(msa):
    for info in msa.generate_debug_log_info():
        log_file_logger.debug(f"{str(info[0] + 1)} {info[1]} {info[2].value} {info[3]}")


def _looks_like_ecomp(path: str) -> bool:
    try:
        with open(path, "rb") as handle:
            header = handle.read(len(HEADER_MAGIC))
    except OSError:
        return False
    return header == HEADER_MAGIC

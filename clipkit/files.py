from enum import Enum
from os import stat
from stat import S_ISREG
from .logger import log_file_logger

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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


_FASTA_BLOCK_BYTES = 1 << 20
_NONPRINTABLE_FASTA_BYTES = bytes(range(32)) + bytes(range(127, 256))
_FASTA_VALIDATION = bytes.maketrans(
    _NONPRINTABLE_FASTA_BYTES, bytes(len(_NONPRINTABLE_FASTA_BYTES))
)


def _iter_fasta_blocks(handle):
    if handle.read(1) != b">":
        raise ValueError("Expected a FASTA header")
    pending = b""
    while block := handle.read(_FASTA_BLOCK_BYTES):
        records = (pending + block).split(b"\n>")
        pending = records.pop()
        yield from records
        if len(pending) > _FASTA_BLOCK_BYTES:
            raise ValueError("FASTA record exceeds the block parser limit")
    yield pending


def _iter_plain_fasta_records(handle):
    for record in _iter_fasta_blocks(handle):
        raw_title, _, sequence = record.partition(b"\n")
        sequence = sequence.translate(_FASTA_VALIDATION, b" \n")
        # Defer unusual whitespace and encodings to the installed Bio parser.
        if b"\0" in sequence or b"\r" in raw_title:
            raise ValueError("Unusual FASTA whitespace or encoding")
        title = raw_title.decode("ascii").rstrip()
        identifier = title.split(None, 1)[0] if title else ""
        yield SeqRecord(Seq(sequence), id=identifier, name=identifier, description=title)


def _read_alignment(path, file_format):
    if file_format == FileFormat.fasta:
        info = stat(path)
        if S_ISREG(info.st_mode) and info.st_size >= _FASTA_BLOCK_BYTES:
            try:
                with open(path, "rb") as handle:
                    return MultipleSeqAlignment(_iter_plain_fasta_records(handle))
            except ValueError:
                pass
    with open(path) as handle:
        return AlignIO.read(handle, file_format.value)


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
        return _read_alignment(input_file_name, file_format), file_format
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
                return _read_alignment(input_file_name, fileFormat), fileFormat
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

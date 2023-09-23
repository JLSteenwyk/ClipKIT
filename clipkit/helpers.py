import re

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from tqdm import tqdm

from .msa import MSA
from .modes import TrimmingMode
from .settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS
from .files import FileFormat
from .stats import TrimmingStats

from enum import Enum


class SeqType(Enum):
    aa = "aa"
    nt = "nt"


def remove_gaps(seq: str, gap_chars: list[str] = DEFAULT_AA_GAP_CHARS) -> str:
    pattern = "|".join([re.escape(char) for char in gap_chars])
    return re.sub(pattern, "", seq)


def get_seq_type(alignment: MultipleSeqAlignment) -> SeqType:
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


def get_alignment_column(alignment: MultipleSeqAlignment, index: int) -> str:
    alignment_column = ""
    alignment_column += alignment[:, index]
    return alignment_column.upper()


def get_gap_chars(seq_type: SeqType) -> list[str]:
    if seq_type == SeqType.nt:
        return DEFAULT_NT_GAP_CHARS
    else:
        return DEFAULT_AA_GAP_CHARS


def report_column_features(
    alignment: MultipleSeqAlignment, index: int, gap_chars: list
) -> tuple[str, float]:
    """
    Count the occurence of each character at a given position
    in an alignment. This information is used to determine
    if the alignment is parsimony informative or not. When
    counting characters, gaps are excluded.
    """
    alignment_column = get_alignment_column(alignment, index)

    column_length = len(alignment_column)

    number_of_gaps = sum([alignment_column.count(char) for char in gap_chars])
    gappyness = number_of_gaps / column_length

    return alignment_column, gappyness


def create_msa(alignment: MultipleSeqAlignment, gap_chars: list[str]=None) -> MSA:
    """
    Create MSA class
    """
    return MSA.from_bio_msa(alignment, gap_chars)


def write_keep_msa(
    msa: MSA, out_file_name: str, out_file_format: FileFormat
) -> None:
    """
    msa is populated with sites that are kept after trimming is finished
    """
    output_msa = msa.to_bio_msa()
    if out_file_format.value == "phylip_relaxed":
        SeqIO.write(output_msa, out_file_name, "phylip-relaxed")
    elif out_file_format.value == "phylip_sequential":
        SeqIO.write(output_msa, out_file_name, "phylip-sequential")
    else:
        SeqIO.write(output_msa, out_file_name, out_file_format.value)


def write_complement(msa: MSA, out_file: str, out_file_format: FileFormat) -> None:
    """
    msa is populated with sites that are trimmed after trimming is finished
    """
    output_msa = msa.complement_to_bio_msa()
    completmentOut = str(out_file) + ".complement"
    if out_file_format.value == "phylip_relaxed":
        SeqIO.write(output_msa, out_file, "phylip-relaxed")
    elif out_file_format.value == "phylip_sequential":
        SeqIO.write(output_msa, out_file, "phylip-sequential")
    SeqIO.write(output_msa, completmentOut, out_file_format.value)


def trim_and_get_stats(
    alignment: MultipleSeqAlignment,
    gaps: float,
    mode: TrimmingMode,
    use_log: bool,
    out_file_name: str,
    complement: bool,
    gap_chars: list,
    quiet: bool,
) -> tuple["TrimRun", "TrimmingStats"]:
    """
    Determines positions to keep or trim and saves these positions on the MSA classes.
    """
    msa = create_msa(alignment)
    msa.trim(mode, gap_threshold=gaps)

    return trim_run, stats
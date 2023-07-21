import re

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from tqdm import tqdm

from .msa import MSA
from .modes import SiteClassificationType, TrimmingMode, trim
from .settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS
from .files import FileFormat
from .write import write_processing_aln, write_output_files_message

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


def count_characters_at_position(alignment_column: str, gap_chars: list) -> dict:
    character_counts = {}

    alignment_column = remove_gaps(alignment_column, gap_chars)
    for char in set(alignment_column):
        character_counts[char] = alignment_column.count(char)
    return character_counts


def determine_site_classification_type(
    character_counts: dict,
) -> SiteClassificationType:
    """
    Determines if a site is parsimony informative or constant.
    A site is parsimony-informative if it contains at least two types of nucleotides
    (or amino acids), and at least two of them occur with a minimum frequency of two.
    https://www.megasoftware.net/web_help_7/rh_parsimony_informative_site.htm

    A site is constant if it contains only one character and that character occurs
    at least twice. https://www.megasoftware.net/web_help_7/rh_constant_site.htm

    A singleton is a site that contains at least two types of characters with, at most,
    one occuring multiple times. https://www.megasoftware.net/web_help_7/rh_singleton_sites.htm
    """
    parsimony_informative_threshold = 2
    counts_gte_threshold = 0

    for count in character_counts.values():
        if count >= 2:
            counts_gte_threshold += 1
        if counts_gte_threshold >= parsimony_informative_threshold:
            return SiteClassificationType.parsimony_informative

    if counts_gte_threshold == 1 and len(character_counts) == 1:
        return SiteClassificationType.constant
    elif counts_gte_threshold == 1 and len(character_counts) > 1:
        return SiteClassificationType.singleton

    return SiteClassificationType.other


def create_keep_and_trim_msas(
    alignment: MultipleSeqAlignment, complement: bool
) -> tuple[MSA, MSA]:
    """
    Creates MSA classes to track sites kept and trimmed
    """

    keep_msa = MSA.from_bio_msa(alignment)
    if complement:
        trim_msa = MSA.from_bio_msa(alignment)
    else:
        trim_msa = None

    return keep_msa, trim_msa


def write_keep_msa(
    keep_msa: MSA, out_file_name: str, out_file_format: FileFormat
) -> None:
    """
    keep_msa is populated with sites that are kept after trimming is finished
    """
    output_msa = keep_msa.to_bio_msa()
    if out_file_format.value == "phylip_relaxed":
        SeqIO.write(output_msa, out_file_name, "phylip-relaxed")
    elif out_file_format.value == "phylip_sequential":
        SeqIO.write(output_msa, out_file_name, "phylip-sequential")
    else:
        SeqIO.write(output_msa, out_file_name, out_file_format.value)


def write_trim_msa(trim_msa: MSA, out_file: str, out_file_format: FileFormat) -> None:
    """
    trim_msa is populated with sites that are trimmed after trimming is finished
    """
    output_msa = trim_msa.to_bio_msa()
    completmentOut = str(out_file) + ".complement"
    if out_file_format.value == "phylip_relaxed":
        SeqIO.write(output_msa, out_file, "phylip-relaxed")
    elif out_file_format.value == "phylip_sequential":
        SeqIO.write(output_msa, out_file, "phylip-sequential")
    SeqIO.write(output_msa, completmentOut, out_file_format.value)


def keep_trim_and_log(
    alignment: MultipleSeqAlignment,
    gaps: float,
    mode: TrimmingMode,
    use_log: bool,
    out_file_name: str,
    complement: bool,
    gap_chars: list,
    quiet: bool,
) -> tuple[MSA, MSA]:
    """
    Determines positions to keep or trim and saves these positions on the MSA classes.
    """
    # initialize MSA classes
    keep_msa, trim_msa = create_keep_and_trim_msas(alignment, complement)

    alignment_length = alignment.get_alignment_length()

    site_classification_counts = dict()
    site_classification_counts[SiteClassificationType.parsimony_informative] = 0
    site_classification_counts[SiteClassificationType.constant] = 0
    site_classification_counts[SiteClassificationType.singleton] = 0
    site_classification_counts[SiteClassificationType.other] = 0

    write_processing_aln()
    for i in tqdm(range(alignment_length), disable=quiet, postfix="trimmer"):
        sequence_at_index, gappyness = report_column_features(alignment, i, gap_chars)

        character_counts = count_characters_at_position(sequence_at_index, gap_chars)

        # determine if a site is parsimony informative, singleton, or constant
        site_classification_type = determine_site_classification_type(character_counts)

        keep_msa, trim_msa = trim(
            gappyness,
            site_classification_type,
            site_classification_counts,
            keep_msa,
            trim_msa,
            i,
            gaps,
            alignment,
            mode,
            use_log,
        )

    # inform user that output files are being written
    write_output_files_message(out_file_name, complement, use_log)
    return keep_msa, trim_msa, site_classification_counts

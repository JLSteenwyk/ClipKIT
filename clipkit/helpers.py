import re

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from typing import Optional

from .ambiguity import is_nucleotide_alphabet
from .msa import MSA
from .modes import AmbiguityHandling, TrimmingMode
from .settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS
from .formats import FileFormat, SeqType
from .stats import TrimmingStats
from .ecomp import write_ecomp

def remove_gaps(seq: str, gap_chars: list[str] = DEFAULT_AA_GAP_CHARS) -> str:
    if isinstance(gap_chars, (list, tuple)) and all(
        isinstance(char, str) and len(char) == 1 for char in gap_chars
    ):
        return seq.translate(str.maketrans("", "", "".join(gap_chars)))
    pattern = "|".join([re.escape(char) for char in gap_chars])
    return re.sub(pattern, "", seq)


def get_seq_type(alignment: MultipleSeqAlignment) -> SeqType:
    seq = str(alignment[0].seq)
    seq = remove_gaps(seq)
    if len(seq) < 200:
        seq = "".join([str(record.seq) for record in alignment])
        seq = remove_gaps(seq)

    sequence_type = SeqType.nt if is_nucleotide_alphabet(set(seq)) else SeqType.aa

    return sequence_type


def get_gap_chars(seq_type: SeqType) -> list[str]:
    if seq_type == SeqType.nt:
        return DEFAULT_NT_GAP_CHARS
    else:
        return DEFAULT_AA_GAP_CHARS


def create_msa(
    alignment: MultipleSeqAlignment,
    gap_chars: list[str] = None,
    threads: int = 1,
    sequence_type: Optional[SeqType] = None,
    ambiguity_handling: AmbiguityHandling = AmbiguityHandling.missing,
) -> MSA:
    """
    Create MSA class
    """
    return MSA.from_bio_msa(
        alignment,
        gap_chars,
        threads,
        sequence_type=sequence_type,
        ambiguity_handling=ambiguity_handling,
    )


def write_msa(
    msa: MSA,
    out_file_name: str,
    out_file_format: FileFormat,
    base_metadata: Optional[dict[str, object]] = None,
) -> None:
    """
    msa is populated with sites that are kept after trimming is finished
    """
    _write_msa_selection(msa, out_file_name, out_file_format, base_metadata)


def write_complement(
    msa: MSA,
    out_file: str,
    out_file_format: FileFormat,
    base_metadata: Optional[dict[str, object]] = None,
) -> None:
    """
    msa is populated with sites that are trimmed after trimming is finished
    """
    _write_msa_selection(
        msa,
        f"{out_file}.complement",
        out_file_format,
        base_metadata,
        complement=True,
    )


def _write_msa_selection(msa, path, file_format, base_metadata, *, complement=False):
    if file_format == FileFormat.fasta and type(msa) is MSA:
        positions = (
            msa._site_positions_to_trim if complement else msa._site_positions_to_keep
        )
        rows = msa._selected_fasta_rows(positions)
        if rows is not None:
            with open(path, "w") as handle:
                for title, residues in rows:
                    if (
                        "\n" in title
                        or "\r" in title
                        or b"\n" in residues
                        or b"\r" in residues
                    ):
                        record = SeqRecord(Seq(residues), id=title, description="")
                        SeqIO.write(record, handle, "fasta")
                    else:
                        _write_fasta_row(handle, title, residues)
            return
    alignment = msa.complement_to_bio_msa() if complement else msa.to_bio_msa()
    _write_alignment(alignment, path, file_format, base_metadata)


def _write_alignment(alignment, path, file_format, base_metadata):
    if file_format == FileFormat.ecomp:
        write_ecomp(alignment, path, base_metadata)
    elif file_format == FileFormat.fasta and alignment.get_alignment_length() >= 4096:
        _write_wrapped_fasta(alignment, path)
    else:
        SeqIO.write(alignment, path, file_format.value.replace("_", "-"))


def _write_wrapped_fasta(alignment, path):
    with open(path, "w") as handle:
        for record in alignment:
            if (
                record.description
                or record.seq is None
                or "\n" in record.id
                or "\r" in record.id
            ):
                SeqIO.write(record, handle, "fasta")
                continue
            residues = bytes(record.seq)
            if not residues.isascii() or b"\n" in residues or b"\r" in residues:
                SeqIO.write(record, handle, "fasta")
                continue
            _write_fasta_row(handle, record.id, residues)


def _write_fasta_row(handle, title, residues):
    # Match Biopython's standard FASTA format, batching the 60-column lines.
    width = 60
    handle.write(f">{title}\n")
    lines, tail = divmod(len(residues), width)
    wrapped = np.empty((lines, width + 1), dtype=np.uint8)
    wrapped[:, :width] = np.frombuffer(
        residues, count=lines * width, dtype=np.uint8
    ).reshape(lines, width)
    wrapped[:, width] = ord("\n")
    handle.write(str(wrapped, "ascii"))
    if tail:
        handle.write(residues[-tail:].decode("ascii") + "\n")

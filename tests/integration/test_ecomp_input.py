import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.ecomp import read_ecomp
from clipkit.modes import TrimmingMode
from clipkit.settings import DEFAULT_AA_GAP_CHARS

here = Path(__file__)
SAMPLES = here.parent / "samples"
EXPECTED = here.parent / "expected"


@pytest.mark.integration
class TestEcompInputs(object):
    def test_kpic_smart_gap_matches_fasta_reference(self):
        input_file = SAMPLES / "12_YIL115C_Anc_2.253_aa_aln.ecomp"
        output_file = (
            "output/12_YIL115C_Anc_2.253_aa_aln.ecomp.clipkit_kpic_smart_gaps"
        )
        expected_file = (
            EXPECTED / "12_YIL115C_Anc_2.253_aa_aln.clipkit_kpic_smart_gaps"
        )

        kwargs = dict(
            input_file=str(input_file),
            output_file=output_file,
            input_file_format="ecomp",
            output_file_format=None,
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9167,
            mode=TrimmingMode.kpic_smart_gap,
            use_log=False,
            gap_characters=DEFAULT_AA_GAP_CHARS,
            quiet=True,
            ends_only=False,
        )
        execute(**kwargs)

        alignment, metadata = read_ecomp(output_file)

        expected_sequences = _parse_fasta(expected_file)
        output_sequences = {record.id: str(record.seq) for record in alignment}

        assert output_sequences == expected_sequences
        assert metadata.get("fallback", {}).get("type") == "gzip"

    def test_explicit_ecomp_output_falls_back_to_fasta(self):
        input_file = SAMPLES / "12_YIL115C_Anc_2.253_aa_aln.ecomp"
        output_file = "output/12_YIL115C_Anc_2.253_aa_aln.ecomp.clipkit_fallback"
        expected_file = (
            EXPECTED / "12_YIL115C_Anc_2.253_aa_aln.clipkit_kpic_smart_gaps"
        )

        kwargs = dict(
            input_file=str(input_file),
            output_file=output_file,
            input_file_format="ecomp",
            output_file_format="ecomp",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9167,
            mode=TrimmingMode.kpic_smart_gap,
            use_log=False,
            gap_characters=DEFAULT_AA_GAP_CHARS,
            quiet=True,
            ends_only=False,
        )
        execute(**kwargs)

        alignment, metadata = read_ecomp(output_file)

        expected_sequences = _parse_fasta(expected_file)
        output_sequences = {record.id: str(record.seq) for record in alignment}

        assert output_sequences == expected_sequences
        assert metadata.get("fallback", {}).get("type") == "gzip"


def _parse_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []
    with path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(current_seq)
    return sequences

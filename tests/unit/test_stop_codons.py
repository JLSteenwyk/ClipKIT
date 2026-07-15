import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from clipkit.modes import StopCodonMode
from clipkit.msa import MSA
from clipkit.settings import DEFAULT_NT_GAP_CHARS


def make_msa(*sequences: str) -> MSA:
    alignment = MultipleSeqAlignment(
        [
            SeqRecord(Seq(sequence), id=f"sequence_{index}")
            for index, sequence in enumerate(sequences)
        ]
    )
    return MSA.from_bio_msa(alignment, DEFAULT_NT_GAP_CHARS)


@pytest.mark.parametrize(
    "mode, expected, terminal_masked, internal_masked",
    [
        (StopCodonMode.terminal, "ATGTAACCT------", 1, 0),
        (StopCodonMode.internal, "ATG---CCTTGA---", 0, 1),
        (StopCodonMode.all, "ATG---CCT------", 1, 1),
    ],
)
def test_masks_requested_stop_codon_types(
    mode, expected, terminal_masked, internal_masked
):
    msa = make_msa("ATGTAACCTTGA---")

    stats = msa.mask_stop_codons(mode)

    assert str(msa.to_bio_msa()[0].seq) == expected
    assert stats.mode is mode
    assert stats.terminal_masked == terminal_masked
    assert stats.internal_masked == internal_masked
    assert stats.total_masked == terminal_masked + internal_masked


def test_handles_lowercase_dna_and_rna_with_trailing_gap_codons():
    msa = make_msa("atgtaaacctga---", "AUGUAAACCUGA---")

    stats = msa.mask_stop_codons(StopCodonMode.all)

    assert [str(record.seq) for record in msa.to_bio_msa()] == [
        "atg---acc------",
        "AUG---ACC------",
    ]
    assert stats.summary == {
        "mode": "all",
        "terminal_masked": 2,
        "internal_masked": 2,
        "total_masked": 4,
    }


def test_ignores_gapped_codons_and_sequences_without_stops():
    msa = make_msa("ATGT-AACCTGA---", "ATGCAAACCCTG---")

    stats = msa.mask_stop_codons(StopCodonMode.all)

    assert [str(record.seq) for record in msa.to_bio_msa()] == [
        "ATGT-AACC------",
        "ATGCAAACCCTG---",
    ]
    assert stats.terminal_masked == 1
    assert stats.internal_masked == 0


def test_terminal_is_the_final_complete_non_gap_codon():
    msa = make_msa("ATGTAAT-A---", "ATGTAA---T-A")

    terminal_stats = msa.mask_stop_codons(StopCodonMode.terminal)

    assert [str(record.seq) for record in msa.to_bio_msa()] == [
        "ATG---T-A---",
        "ATG------T-A",
    ]
    assert terminal_stats.terminal_masked == 2


def test_masking_invalidates_cached_gap_statistics():
    msa = make_msa("ATGTAA", "ATGCAA")
    np.testing.assert_equal(msa.site_gappyness, np.zeros(6))

    msa.mask_stop_codons(StopCodonMode.terminal)

    np.testing.assert_equal(msa.site_gappyness, np.array([0, 0, 0, 0.5, 0.5, 0.5]))


def test_mask_character_becomes_an_effective_gap_when_not_configured():
    alignment = MultipleSeqAlignment(
        [
            SeqRecord(Seq("ATGTAA---"), id="stop"),
            SeqRecord(Seq("ATGCAA---"), id="control"),
        ]
    )
    msa = MSA.from_bio_msa(alignment, gap_chars=["?"])

    stats = msa.mask_stop_codons(StopCodonMode.terminal)

    assert stats.terminal_masked == 1
    assert "-" in msa.gap_chars
    np.testing.assert_equal(msa.site_gappyness[3:], np.array([0.5, 0.5, 0.5, 1, 1, 1]))


def test_requires_complete_codon_columns():
    msa = make_msa("ATGTA")

    with pytest.raises(ValueError, match="alignment length divisible by 3"):
        msa.mask_stop_codons(StopCodonMode.all)

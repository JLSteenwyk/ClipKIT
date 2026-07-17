import numpy as np
import pytest
import random
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


def test_terminal_stops_can_precede_different_numbers_of_trailing_gap_codons():
    msa = make_msa("ATGTAA------", "ATGCAATAG---")

    terminal_stats = msa.mask_stop_codons(StopCodonMode.terminal)

    assert [str(record.seq) for record in msa.to_bio_msa()] == [
        "ATG---------",
        "ATGCAA------",
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


def _scalar_stop_codon_mask(sequences, mode):
    gap_chars = {char.upper() for char in DEFAULT_NT_GAP_CHARS}
    expected = [list(sequence) for sequence in sequences]
    terminal_masked = 0
    internal_masked = 0

    for row_index, sequence in enumerate(sequences):
        codons = [sequence[index : index + 3] for index in range(0, len(sequence), 3)]
        complete_positions = [
            index
            for index, codon in enumerate(codons)
            if not any(char.upper() in gap_chars for char in codon)
        ]
        terminal_position = complete_positions[-1] if complete_positions else None

        for codon_index, codon in enumerate(codons):
            is_complete = codon_index in complete_positions
            is_stop = is_complete and codon.upper() in {
                "TAA",
                "TAG",
                "TGA",
                "UAA",
                "UAG",
                "UGA",
            }
            is_terminal = codon_index == terminal_position
            selected = is_stop and (
                mode is StopCodonMode.all
                or (mode is StopCodonMode.terminal and is_terminal)
                or (mode is StopCodonMode.internal and not is_terminal)
            )
            if not selected:
                continue

            start = codon_index * 3
            expected[row_index][start : start + 3] = "---"
            if is_terminal:
                terminal_masked += 1
            else:
                internal_masked += 1

    return (
        ["".join(sequence) for sequence in expected],
        terminal_masked,
        internal_masked,
    )


@pytest.mark.parametrize("mode", list(StopCodonMode))
@pytest.mark.parametrize("seed", range(8))
def test_stop_codon_masking_matches_randomized_scalar_reference(mode, seed):
    rng = random.Random(seed)
    choices = ["ATG", "CAA", "TTC", "GGA", "TAA", "TAG", "TGA", "---", "A-G"]
    sequences = []
    for row in range(18):
        codons = [rng.choice(choices) for _ in range(40)]
        if row % 3 == 0:
            codons = [codon.lower() for codon in codons]
        if row % 4 == 0:
            codons[-3:] = ["---", "---", "---"]
        sequences.append("".join(codons))

    expected, terminal_masked, internal_masked = _scalar_stop_codon_mask(
        sequences, mode
    )
    msa = make_msa(*sequences)

    stats = msa.mask_stop_codons(mode)

    assert [str(record.seq) for record in msa.to_bio_msa()] == expected
    assert stats.terminal_masked == terminal_masked
    assert stats.internal_masked == internal_masked

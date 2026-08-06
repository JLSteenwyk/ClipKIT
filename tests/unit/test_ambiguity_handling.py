import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from clipkit.ambiguity import ambiguity_symbols
from clipkit.helpers import SeqType, get_seq_type
from clipkit.modes import AmbiguityHandling, TrimmingMode
from clipkit.msa import MSA
from clipkit.site_classification import SiteClassificationType


def _msa(column, *, sequence_type="nt", handling="missing", gap_chars=None):
    records = np.array([[character] for character in column], dtype="U1")
    return MSA(
        [{"id": str(index), "description": str(index)} for index in range(len(column))],
        records,
        gap_chars=gap_chars or ["-"],
        sequence_type=sequence_type,
        ambiguity_handling=handling,
    )


@pytest.mark.parametrize(
    "handling, expected_entropy, expected_gappyness",
    [
        (AmbiguityHandling.missing, 0.9183, 0.25),
        (AmbiguityHandling.fractional, 0.9544, 0.0),
        (AmbiguityHandling.literal, 0.9464, 0.0),
    ],
)
def test_nucleotide_ambiguity_policies(handling, expected_entropy, expected_gappyness):
    msa = _msa("AAGR", handling=handling)

    assert msa.site_entropy[0] == expected_entropy
    assert msa.site_gappyness[0] == expected_gappyness
    assert msa.site_gap_fraction[0] == 0.0
    assert msa.site_ambiguity[0] == 0.25
    assert msa.site_resolved_fraction[0] == 0.75


@pytest.mark.parametrize(
    "handling, expected",
    [
        (AmbiguityHandling.missing, SiteClassificationType.constant),
        (AmbiguityHandling.fractional, SiteClassificationType.constant),
        (AmbiguityHandling.literal, SiteClassificationType.parsimony_informative),
    ],
)
def test_ambiguity_never_creates_pi_except_in_literal_mode(handling, expected):
    msa = _msa("AARR", handling=handling)

    assert msa.site_classification_types[0] == expected


def test_fractional_protein_and_rna_expansions():
    protein = _msa("BBDN", sequence_type="aa", handling="fractional")
    rna = _msa("UCYY", sequence_type="nt", handling="fractional")

    assert protein.column_character_frequencies == [{"D": 2.0, "N": 2.0}]
    assert rna.column_character_frequencies == [{"C": 2.0, "U": 2.0}]
    assert protein.site_entropy[0] == 1.0
    assert rna.site_entropy[0] == 1.0


def test_configured_gap_characters_take_precedence_over_fractional_expansion():
    msa = _msa("AANN", handling="fractional", gap_chars=["-", "N"])

    assert msa.column_character_frequencies == [{"A": 2}]
    assert msa.site_gap_fraction[0] == 0.5
    assert msa.site_gappyness[0] == 0.5


@pytest.mark.parametrize("mode", [TrimmingMode.entropy, TrimmingMode.composition_bias])
def test_all_missing_site_is_trimmed_by_state_scoring_modes(mode):
    msa = _msa("RRRR", handling="missing")

    msa.trim(mode=mode, gap_threshold=0.8)

    np.testing.assert_equal(msa._site_positions_to_trim, np.array([0]))
    assert msa.site_analyzable_fraction[0] == 0.0


def test_all_missing_site_is_trimmed_by_heterotachy_mode():
    msa = _msa("RRRR", handling="missing")

    msa.trim(mode=TrimmingMode.heterotachy, gap_threshold=0.8, guide_tree=Tree())

    np.testing.assert_equal(msa._site_positions_to_trim, np.array([0]))


def test_ambiguity_handling_does_not_rewrite_alignment_characters():
    alignment = MultipleSeqAlignment(
        [
            SeqRecord(Seq(character), id=str(index))
            for index, character in enumerate("RrYN")
        ]
    )
    msa = MSA.from_bio_msa(
        alignment,
        gap_chars=["-"],
        sequence_type="nt",
        ambiguity_handling="fractional",
    )

    assert [str(record.seq) for record in msa.to_bio_msa()] == ["R", "r", "Y", "N"]
    assert msa.site_ambiguity[0] == 1.0
    assert msa.column_character_frequencies == [
        {"A": 1.25, "C": 0.75, "G": 1.25, "T": 0.75}
    ]


def test_iupac_rich_alignment_is_auto_detected_as_nucleotide():
    alignment = MultipleSeqAlignment(
        [
            SeqRecord(Seq("ACGTRYSWKMBDHVN"), id="one"),
            SeqRecord(Seq("NNNNACGTURY----"), id="two"),
        ]
    )

    assert get_seq_type(alignment) == SeqType.nt

    msa = MSA.from_bio_msa(alignment, gap_chars=["-"])
    assert msa.sequence_type == "nt"
    assert msa.site_ambiguity[4] == 0.5


def test_amino_acid_alignment_is_not_misdetected_as_nucleotide():
    alignment = MultipleSeqAlignment(
        [SeqRecord(Seq("MPEPTIDE"), id="one"), SeqRecord(Seq("MPEPTVDE"), id="two")]
    )

    assert get_seq_type(alignment) == SeqType.aa


def test_supported_iupac_ambiguity_symbols_are_sequence_type_specific():
    assert ambiguity_symbols("nt") == frozenset("RYSWKMBDHVNX")
    assert ambiguity_symbols("aa") == frozenset("BZJX")

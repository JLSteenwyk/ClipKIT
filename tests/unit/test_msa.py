import pytest
import numpy as np
import math
import random

from Bio import AlignIO
from clipkit.guide_tree import build_parsimony_guide_tree
from clipkit.msa import MSA, _column_character_counts
from clipkit.modes import TrimmingMode
from clipkit.site_classification import determine_site_classification_type


def get_biopython_msa(file_path, file_format="fasta"):
    return AlignIO.read(open(file_path), file_format)


class TestMSA(object):
    def test_clipkit_msa_from_bio_msa(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        assert msa.header_info == [
            {"id": "1", "name": "1", "description": "1"},
            {"id": "2", "name": "2", "description": "2"},
            {"id": "3", "name": "3", "description": "3"},
            {"id": "4", "name": "4", "description": "4"},
            {"id": "5", "name": "5", "description": "5"},
        ]
        expected_seq_records = np.array(
            [
                ["A", "-", "G", "T", "A", "T"],
                ["A", "-", "G", "-", "A", "T"],
                ["A", "-", "G", "-", "T", "A"],
                ["A", "G", "A", "-", "T", "A"],
                ["A", "C", "a", "-", "T", "-"],
            ]
        )
        np.testing.assert_equal(msa.seq_records, expected_seq_records)

    def test_to_bio_msa_preserves_sequences_and_descriptions(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)

        converted = msa.to_bio_msa()

        assert [str(record.seq) for record in converted] == [
            str(record.seq) for record in bio_msa
        ]
        assert [record.id for record in converted] == [
            record.description for record in bio_msa
        ]

    def test_to_bio_msa_handles_empty_sequence_rows(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        msa.trim(site_positions_to_trim=np.arange(msa.original_length))

        converted = msa.to_bio_msa()

        assert [str(record.seq) for record in converted] == [""] * len(bio_msa)

    def test_trim_by_provided_site_positions_np_array(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        sites_to_trim = np.array([1, 4])
        msa.trim(site_positions_to_trim=sites_to_trim)
        expected_sites_kept = np.array(
            [
                ["A", "G", "T", "T"],
                ["A", "G", "-", "T"],
                ["A", "G", "-", "A"],
                ["A", "A", "-", "A"],
                ["A", "a", "-", "-"],
            ]
        )
        np.testing.assert_equal(msa.sites_kept, expected_sites_kept)

    def test_trim_by_provided_site_positions_list(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        sites_to_trim = [1, 4]
        msa.trim(site_positions_to_trim=sites_to_trim)
        expected_sites_kept = np.array(
            [
                ["A", "G", "T", "T"],
                ["A", "G", "-", "T"],
                ["A", "G", "-", "A"],
                ["A", "A", "-", "A"],
                ["A", "a", "-", "-"],
            ]
        )
        np.testing.assert_equal(msa.sites_kept, expected_sites_kept)

    @pytest.mark.parametrize(
        "sites_to_trim, expected",
        [
            (
                [0],
                np.array(
                    [
                        ["T", "A", "T"],
                        ["-", "A", "T"],
                        ["-", "T", "A"],
                        ["-", "T", "A"],
                        ["-", "T", "-"],
                    ]
                ),
            ),
            (
                [2],
                np.array(
                    [
                        ["T", "A", "T"],
                        ["-", "A", "T"],
                        ["-", "T", "A"],
                        ["-", "T", "A"],
                        ["-", "T", "-"],
                    ]
                ),
            ),
            (
                [3],
                np.array(
                    [
                        ["A", "-", "G"],
                        ["A", "-", "G"],
                        ["A", "-", "G"],
                        ["A", "G", "A"],
                        ["A", "C", "a"],
                    ]
                ),
            ),
            (
                [5],
                np.array(
                    [
                        ["A", "-", "G"],
                        ["A", "-", "G"],
                        ["A", "-", "G"],
                        ["A", "G", "A"],
                        ["A", "C", "a"],
                    ]
                ),
            ),
            (
                [0, 1, 2, 3, 4, 5],
                np.array(
                    [
                        [],
                        [],
                        [],
                        [],
                        [],
                    ],
                    dtype=object,
                ),
            ),
        ],
    )
    def test_trim_codons(self, sites_to_trim, expected):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        msa.trim(site_positions_to_trim=sites_to_trim, codon=True)
        np.testing.assert_equal(msa.trimmed, expected)

    def test_entropy_mode_trims_high_entropy_sites(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        msa.trim(mode=TrimmingMode.entropy, gap_threshold=0.95)
        np.testing.assert_equal(msa._site_positions_to_trim, np.array([1, 2, 4, 5]))

    def test_gappyout_threshold_on_simple_alignment(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        assert msa.determine_gappyout_gap_threshold() == 0.4

    def test_gappyout_mode_trims_high_gap_outliers(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        msa.trim(mode=TrimmingMode.gappyout, gap_threshold=0.4)

        expected_sites_kept = np.array(
            [
                ["A", "G", "A", "T"],
                ["A", "G", "A", "T"],
                ["A", "G", "T", "A"],
                ["A", "A", "T", "A"],
                ["A", "a", "T", "-"],
            ]
        )
        np.testing.assert_equal(msa.sites_kept, expected_sites_kept)

    def test_block_gappy_mode_trims_only_contiguous_high_gap_sites(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        msa.trim(mode=TrimmingMode.block_gappy, gap_threshold=0.6)
        np.testing.assert_equal(msa._site_positions_to_trim, np.array([], dtype=int))

    def test_composition_bias_mode_trims_strongly_dominated_sites(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        msa.trim(mode=TrimmingMode.composition_bias, gap_threshold=0.9)
        np.testing.assert_equal(msa._site_positions_to_trim, np.array([0, 3]))

    def test_heterotachy_mode_trims_sites_with_high_clade_entropy_variation(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        guide_tree = build_parsimony_guide_tree(bio_msa)
        msa.trim(
            mode=TrimmingMode.heterotachy,
            gap_threshold=0.8,
            guide_tree=guide_tree,
        )
        np.testing.assert_equal(msa._site_positions_to_trim, np.array([1, 2]))


@pytest.mark.parametrize(
    "seq_records",
    [
        np.array([list("Aa-?X"), list("AA--x"), list("C?-ax")], dtype="U1"),
        np.array(
            [["A", "Ā", "😀"], ["Ω", "A", "😀"], ["Ā", "Ω", "A"]],
            dtype="U1",
        ),
    ],
)
def test_column_character_counts_match_per_column_unique(seq_records):
    states, counts = _column_character_counts(seq_records)
    expected_states = np.unique(seq_records)
    expected_counts = np.array(
        [np.count_nonzero(seq_records == state, axis=0) for state in expected_states]
    )

    np.testing.assert_equal(states, expected_states)
    np.testing.assert_equal(counts, expected_counts)


def test_count_backed_properties_match_per_column_reference():
    seq_records = np.array(
        [
            list("Aa-?XZ"),
            list("AA--xZ"),
            list("C?-axZ"),
            list("CA?aXZ"),
        ],
        dtype="U1",
    )
    gap_chars = ["-", "?", "X", "x"]
    msa = MSA(
        [{"id": str(idx)} for idx in range(len(seq_records))],
        seq_records,
        gap_chars=gap_chars,
        requires_uppercase_normalization=True,
    )

    normalized = np.char.upper(seq_records)
    expected_frequencies = []
    expected_entropy = []
    expected_bias = []
    gap_chars_upper = {gap.upper() for gap in gap_chars}
    for column in normalized.T:
        states, counts = np.unique(column, return_counts=True)
        frequencies = dict(zip(states, counts))
        for gap_char in gap_chars:
            frequencies.pop(gap_char, None)
        expected_frequencies.append(frequencies)

        entropy_counts = [
            count
            for state, count in zip(states, counts)
            if state.upper() not in gap_chars_upper
        ]
        if len(entropy_counts) <= 1:
            expected_entropy.append(0.0)
        else:
            total = float(sum(entropy_counts))
            probabilities = [count / total for count in entropy_counts]
            raw_entropy = -sum(
                probability * math.log2(probability)
                for probability in probabilities
                if probability > 0.0
            )
            expected_entropy.append(raw_entropy / math.log2(len(entropy_counts)))

        bias_counts = np.array(list(frequencies.values()), dtype=float)
        total = bias_counts.sum()
        if total == 0:
            expected_bias.append(0.0)
        elif len(bias_counts) == 1:
            expected_bias.append(1.0)
        else:
            probabilities = bias_counts / total
            dominance = float(np.sum(np.square(probabilities)))
            min_dominance = 1.0 / float(len(bias_counts))
            expected_bias.append((dominance - min_dominance) / (1.0 - min_dominance))

    assert msa.column_character_frequencies == expected_frequencies
    np.testing.assert_equal(
        msa.site_classification_types,
        np.array(
            [
                determine_site_classification_type(frequencies)
                for frequencies in expected_frequencies
            ]
        ),
    )
    np.testing.assert_equal(
        msa.site_gappyness,
        np.around(np.isin(seq_records, gap_chars).mean(axis=0), decimals=4),
    )
    np.testing.assert_equal(msa.site_entropy, np.around(expected_entropy, decimals=4))
    np.testing.assert_equal(
        msa.site_composition_bias, np.around(expected_bias, decimals=4)
    )


@pytest.mark.parametrize("seed", range(40))
def test_entropy_and_composition_bias_match_randomized_scalar_reference(seed):
    rng = random.Random(seed)
    row_count = rng.randint(1, 40)
    column_count = rng.randint(1, 160)
    alphabet = "ACDEFGHIKLMNPQRSTVWYacgt-?*Xx"
    seq_records = np.array(
        [[rng.choice(alphabet) for _ in range(column_count)] for _ in range(row_count)],
        dtype="U1",
    )
    gap_chars = ["-", "?", "*", "X", "x"]
    msa = MSA(
        [{"id": str(index)} for index in range(row_count)],
        seq_records,
        gap_chars=gap_chars,
        requires_uppercase_normalization=bool(
            np.any(seq_records != np.char.upper(seq_records))
        ),
    )

    normalized = np.char.upper(seq_records)
    entropy_gap_chars = {char.upper() for char in gap_chars}
    expected_entropy = []
    expected_bias = []
    for column in normalized.T:
        states, counts = np.unique(column, return_counts=True)

        entropy_counts = [
            count
            for state, count in zip(states, counts)
            if state not in entropy_gap_chars
        ]
        if len(entropy_counts) <= 1:
            expected_entropy.append(0.0)
        else:
            total = float(sum(entropy_counts))
            probabilities = [count / total for count in entropy_counts]
            raw_entropy = -sum(
                probability * math.log2(probability)
                for probability in probabilities
                if probability > 0.0
            )
            expected_entropy.append(raw_entropy / math.log2(len(entropy_counts)))

        bias_counts = np.array(
            [count for state, count in zip(states, counts) if state not in gap_chars],
            dtype=float,
        )
        total = bias_counts.sum()
        if total == 0:
            expected_bias.append(0.0)
        elif len(bias_counts) == 1:
            expected_bias.append(1.0)
        else:
            probabilities = bias_counts / total
            dominance = float(np.sum(np.square(probabilities)))
            min_dominance = 1.0 / float(len(bias_counts))
            expected_bias.append((dominance - min_dominance) / (1.0 - min_dominance))

    np.testing.assert_equal(msa.site_entropy, np.around(expected_entropy, decimals=4))
    np.testing.assert_equal(
        msa.site_composition_bias, np.around(expected_bias, decimals=4)
    )

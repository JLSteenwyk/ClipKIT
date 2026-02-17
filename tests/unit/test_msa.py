import pytest
import numpy as np

from Bio import AlignIO
from clipkit.msa import MSA
from clipkit.modes import TrimmingMode


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

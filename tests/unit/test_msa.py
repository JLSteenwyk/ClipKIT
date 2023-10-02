import pytest
import numpy as np

from Bio import AlignIO
from clipkit.msa import MSA

def get_biopython_msa(file_path, file_format="fasta"):
    return AlignIO.read(open(file_path), file_format)


class TestMSA(object):
    def test_clipkit_msa_from_bio_msa(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        assert msa.header_info == [{'id': '1', 'name': '1', 'description': '1'}, {'id': '2', 'name': '2', 'description': '2'}, {'id': '3', 'name': '3', 'description': '3'}, {'id': '4', 'name': '4', 'description': '4'}, {'id': '5', 'name': '5', 'description': '5'}]
        expected_seq_records = np.array([
            ['A', '-', 'G', 'T', 'A', 'T'],
            ['A', '-', 'G', '-', 'A', 'T'],
            ['A', '-', 'G', '-', 'T', 'A'],
            ['A', 'G', 'A', '-', 'T', 'A'],
            ['A', 'C', 'a', '-', 'T', '-']
        ])
        np.testing.assert_equal(msa.seq_records, expected_seq_records)

    def test_trim_by_provided_site_positions_np_array(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        sites_to_trim = np.array([1, 4])
        msa.trim(site_positions_to_trim=sites_to_trim)
        expected_sites_kept = np.array([
            ['A', 'G', 'T', 'T'],
            ['A', 'G', '-', 'T'],
            ['A', 'G', '-', 'A'],
            ['A', 'A', '-', 'A'],
            ['A', 'a', '-', '-']
        ])
        np.testing.assert_equal(msa.sites_kept, expected_sites_kept)

    def test_trim_by_provided_site_positions_list(self):
        bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
        msa = MSA.from_bio_msa(bio_msa)
        sites_to_trim = [1, 4]
        msa.trim(site_positions_to_trim=sites_to_trim)
        expected_sites_kept = np.array([
            ['A', 'G', 'T', 'T'],
            ['A', 'G', '-', 'T'],
            ['A', 'G', '-', 'A'],
            ['A', 'A', '-', 'A'],
            ['A', 'a', '-', '-']
        ])
        np.testing.assert_equal(msa.sites_kept, expected_sites_kept)

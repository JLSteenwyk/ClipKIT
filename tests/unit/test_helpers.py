import pytest
import pytest_mock
from pathlib import Path


import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from clipkit.helpers import count_characters_at_position
from clipkit.helpers import report_column_features
from clipkit.helpers import determine_site_classification_type
from clipkit.helpers import create_keep_and_trim_msas
from clipkit.helpers import write_trim_msa
from clipkit.helpers import write_keep_msa
from clipkit.helpers import SeqType
from clipkit.files import FileFormat
from clipkit.modes import SiteClassificationType, TrimmingMode, trim, should_keep_site
from clipkit.settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS


here = Path(__file__)


@pytest.fixture
def sample_msa():
    return MultipleSeqAlignment(
        [
            SeqRecord(
                seq=Seq("['A']"),
                id="1",
                name="<unknown name>",
                description="",
                dbxrefs=[],
            ),
            SeqRecord(
                seq=Seq("['A']"),
                id="2",
                name="<unknown name>",
                description="",
                dbxrefs=[],
            ),
            SeqRecord(
                seq=Seq("['A']"),
                id="3",
                name="<unknown name>",
                description="",
                dbxrefs=[],
            ),
            SeqRecord(
                seq=Seq("['A']"),
                id="4",
                name="<unknown name>",
                description="",
                dbxrefs=[],
            ),
            SeqRecord(
                seq=Seq("['A']"),
                id="5",
                name="<unknown name>",
                description="",
                dbxrefs=[],
            ),
        ]
    )


class TestCountCharactersAtPosition(object):
    def test_gives_count_for_each_char(self):
        ## setup
        s = "ACTTTGGG"

        ## execution
        res = count_characters_at_position(s, DEFAULT_NT_GAP_CHARS)

        ## check results
        # test that each character has an associated key
        for char in s:
            assert char in res.keys()

        # test that the len of the res is equal to the
        # number of unique string characters
        assert len(res) == len(set(s))


class TestGetSequenceAtPositionAndReportFeatures(object):
    def test_gets_sequence_and_gappyness(self):
        ## setup
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        i = 5

        ## execution
        seq, gappyness = report_column_features(alignment, i, DEFAULT_AA_GAP_CHARS)

        ## check results
        # test output types
        assert isinstance(seq, str)
        assert isinstance(gappyness, float)


class TestParsimonyInformativeOrConstant(object):
    def test_determine_site_classification_type(self):
        ## set up
        num_occurences_pi = {"A": 5, "T": 10, "G": 2, "C": 4}
        num_occurences_const = {"T": 10}
        num_occurences_singleton = {"A": 1, "T": 2}
        num_occurences_other = {"A": 1}

        ## execution
        res_num_occurences_pi = determine_site_classification_type(num_occurences_pi)
        res_num_occurences_const = determine_site_classification_type(
            num_occurences_const
        )
        res_num_occurences_singleton = determine_site_classification_type(
            num_occurences_singleton
        )
        res_num_occurences_other = determine_site_classification_type(
            num_occurences_other
        )

        ## check results
        assert res_num_occurences_pi is SiteClassificationType.parsimony_informative
        assert res_num_occurences_const is SiteClassificationType.constant
        assert res_num_occurences_singleton is SiteClassificationType.singleton
        assert res_num_occurences_other is SiteClassificationType.other


class TestPopulateEmptyKeepDAndTrimD(object):
    def test_create_keep_and_trim_msas(self):
        ## set up
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")

        ## execution
        keep_msa, trim_msa = create_keep_and_trim_msas(alignment, True)

        ## check results
        expected_keep_data = {
            "1": np.array([b"", b"", b"", b"", b"", b""]),
            "2": np.array([b"", b"", b"", b"", b"", b""]),
            "3": np.array([b"", b"", b"", b"", b"", b""]),
            "4": np.array([b"", b"", b"", b"", b"", b""]),
            "5": np.array([b"", b"", b"", b"", b"", b""]),
        }
        expected_trim_data = {
            "1": np.array([b"", b"", b"", b"", b"", b""]),
            "2": np.array([b"", b"", b"", b"", b"", b""]),
            "3": np.array([b"", b"", b"", b"", b"", b""]),
            "4": np.array([b"", b"", b"", b"", b"", b""]),
            "5": np.array([b"", b"", b"", b"", b"", b""]),
        }

        assert expected_keep_data.keys() == keep_msa._data.keys()
        assert all(
            np.array_equal(expected_keep_data[key], keep_msa._data[key])
            for key in expected_keep_data
        )
        assert expected_trim_data.keys() == trim_msa._data.keys()
        assert all(
            np.array_equal(expected_trim_data[key], trim_msa._data[key])
            for key in expected_trim_data
        )


class TestWriteKeepD(object):
    def test_write_keep_msa_writes_file(self, mocker, sample_msa):
        ## set up
        alignment = AlignIO.read(f"{here.parent}/examples/single_site.fa", "fasta")
        keep_msa, _ = create_keep_and_trim_msas(alignment, True)

        out_file = "output_file_name.fa"
        out_file_format = FileFormat.fasta
        mock_msa = mocker.patch("clipkit.msa.MultipleSeqAlignment")
        mock_msa.return_value = sample_msa
        mock_write = mocker.patch("clipkit.helpers.SeqIO.write")

        ## execution
        write_keep_msa(keep_msa, out_file, out_file_format)

        ## check results
        mock_write.assert_called_once_with(sample_msa, out_file, out_file_format.value)


class TestWriteTrimD(object):
    def test_write_trim_msa_calls_seqio_write(self, mocker, sample_msa):
        ## set up
        alignment = AlignIO.read(f"{here.parent}/examples/single_site.fa", "fasta")
        _, trim_msa = create_keep_and_trim_msas(alignment, True)

        out_file = "output_file_name.fa"
        out_file_format = FileFormat.fasta
        mock_msa = mocker.patch("clipkit.msa.MultipleSeqAlignment")
        mock_msa.return_value = sample_msa
        mock_write = mocker.patch("clipkit.helpers.SeqIO.write")

        ## execution
        write_trim_msa(trim_msa, out_file, out_file_format)

        ## check results
        expected_completmentOut = f"{out_file}.complement"
        mock_write.assert_called_once_with(
            sample_msa, expected_completmentOut, out_file_format.value
        )

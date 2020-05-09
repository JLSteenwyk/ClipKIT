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
from clipkit.helpers import get_sequence_at_position_and_report_features
from clipkit.helpers import parsimony_informative_or_constant
from clipkit.helpers import populate_empty_keepD_and_trimD
from clipkit.helpers import join_keepD_and_trimD
from clipkit.helpers import write_trimD
from clipkit.helpers import write_keepD
from clipkit.files import FileFormat

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
        res = count_characters_at_position(s)

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
        i = int(5)

        ## execution
        seq, gappyness = get_sequence_at_position_and_report_features(alignment, i)

        ## check results
        # test output types
        assert isinstance(seq, str)
        assert isinstance(gappyness, float)


class TestParsimonyInformativeOrConstant(object):
    def test_parsimony_informative_or_constant(self):
        ## set up
        # pi = parsimony informative
        num_occurences_pi = {"A": 5, "T": 10, "G": 2, "C": 4}
        # npi = not parsimony informative
        num_occurences_npi = {"A": 1, "T": 10, "G": 1}
        # Const = constant
        num_occurences_const = {"A": 10}
        # nConst = not constant
        num_occurences_nconst = {"A": 1}

        ## execution
        # result is True and False
        (
            is_parsimony_informative,
            constant_site_holder_is_pi,
        ) = parsimony_informative_or_constant(num_occurences_pi)
        # result is False and False
        (
            is_not_parsimony_informative,
            constant_site_holder_is_npi,
        ) = parsimony_informative_or_constant(num_occurences_npi)
        # result is False and True
        is_not_pi_0, is_constant_site = parsimony_informative_or_constant(
            num_occurences_const
        )
        # result is False and False
        is_not_pi_1, is_not_constant_site = parsimony_informative_or_constant(
            num_occurences_nconst
        )

        ## check results
        assert is_parsimony_informative == True and constant_site_holder_is_pi == False
        assert (
            is_not_parsimony_informative == False
            and constant_site_holder_is_npi == False
        )
        assert is_not_pi_0 == False and is_constant_site == True
        assert is_not_pi_1 == False and is_not_constant_site == False


class TestPopulateEmptyKeepDAndTrimD(object):
    def test_populate_empty_keepD_and_trimD(self):
        ## set up
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")

        ## execution
        keepD, trimD = populate_empty_keepD_and_trimD(alignment)

        ## check results
        expected_keepD = {
            "1": np.empty([6], dtype=str),
            "2": np.empty([6], dtype=str),
            "3": np.empty([6], dtype=str),
            "4": np.empty([6], dtype=str),
            "5": np.empty([6], dtype=str),
        }
        expected_trimD = {
            "1": np.empty([6], dtype=str),
            "2": np.empty([6], dtype=str),
            "3": np.empty([6], dtype=str),
            "4": np.empty([6], dtype=str),
            "5": np.empty([6], dtype=str),
        }

        assert expected_keepD.keys() == keepD.keys()
        assert all(
            np.array_equal(expected_keepD[key], keepD[key]) for key in expected_keepD
        )
        assert expected_trimD.keys() == trimD.keys()
        assert all(
            np.array_equal(expected_trimD[key], trimD[key]) for key in expected_trimD
        )


class TestJoinKeepDAndTrimD(object):
    def test_join_keepD_and_trimD(self):
        ## set up
        keepD = {
            "1": ["A", "-", "G", "T", "A", "T"],
            "2": ["A", "-", "G", "-", "A", "T"],
            "3": ["A", "-", "G", "-", "T", "A"],
            "4": ["A", "G", "A", "-", "T", "A"],
            "5": ["A", "C", "a", "-", "T", "-"],
        }
        trimD = {"1": [], "2": [], "3": [], "4": [], "5": []}

        ## execution
        keepD, trimD = join_keepD_and_trimD(keepD, trimD)

        ## check results
        expected_keepD = {
            "1": "A-GTAT",
            "2": "A-G-AT",
            "3": "A-G-TA",
            "4": "AGA-TA",
            "5": "ACa-T-",
        }
        expected_trimD = {"1": "", "2": "", "3": "", "4": "", "5": ""}

        assert expected_keepD == keepD
        assert expected_trimD == trimD


class TestWriteKeepD(object):
    def test_write_keepD_writes_file(self, mocker, sample_msa):
        ## set up
        keepD = {"1": ["A"], "2": ["A"], "3": ["A"], "4": ["A"], "5": ["A"]}
        out_file = "output_file_name.fa"
        out_file_format = FileFormat.fasta
        mock_msa = mocker.patch("clipkit.helpers.MultipleSeqAlignment")
        mock_msa.return_value = sample_msa
        mock_write = mocker.patch("clipkit.helpers.SeqIO.write")

        ## execution
        write_keepD(keepD, out_file, out_file_format)

        ## check results
        mock_write.assert_called_once_with(sample_msa, out_file, out_file_format.value)


class TestWriteTrimD(object):
    def test_write_trimD_calls_seqio_write(self, mocker, sample_msa):
        ## set up
        trimD = {"1": ["A"], "2": ["A"], "3": ["A"], "4": ["A"], "5": ["A"]}
        out_file = "output_file_name.fa"
        out_file_format = FileFormat.fasta
        mock_msa = mocker.patch("clipkit.helpers.MultipleSeqAlignment")
        mock_msa.return_value = sample_msa
        mock_write = mocker.patch("Bio.SeqIO.write")

        ## execution
        write_trimD(trimD, out_file, out_file_format)

        ## check results
        expected_completmentOut = f"{out_file}.complement"
        mock_write.assert_called_once_with(
            sample_msa, expected_completmentOut, out_file_format.value
        )

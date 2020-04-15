import pytest
import pytest_mock

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from pathlib import Path

from clipkit.helpers import count_characters_at_position
from clipkit.helpers import get_sequence_at_position_and_report_features
from clipkit.helpers import determine_if_parsimony_informative
from clipkit.helpers import populate_empty_keepD_and_trimD
from clipkit.helpers import join_keepD_and_trimD
from clipkit.helpers import write_trimD
from clipkit.helpers import write_keepD
from clipkit.files import FileFormat


here = Path(__file__)

class TestCountCharactersAtPosition(object):

    def test_gives_count_for_each_char(self):
        ## setup
        s = 'ACTTTGGG'

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

class TestDetermineIfParsimonyInformative(object):

    def test_determine_if_parsimony_informative(self):
        ## set up
        # pi = parsimony informative
        num_occurences_pi = {'A':5, 'T':10, 'G':2, 'C':4}
        # npi = not parsimony informative
        num_occurences_npi = {'A':1, 'T':10, 'G':1}

        ## execution
        # result is True
        is_parsimony_informative = determine_if_parsimony_informative(num_occurences_pi)
        # result is False
        is_not_parsimony_informative = determine_if_parsimony_informative(num_occurences_npi)

        ## check results
        assert is_parsimony_informative == True
        assert is_not_parsimony_informative == False

class TestPopulateEmptyKeepDAndTrimD(object):

    def test_populate_empty_keepD_and_trimD(self):
        ## set up
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")

        ## execution
        keepD, trimD, logArr = populate_empty_keepD_and_trimD(alignment)

        ## check results
        expected_keepD = {'1':[], '2':[], '3':[], '4':[], '5':[]}
        expected_trimD = {'1':[], '2':[], '3':[], '4':[], '5':[]}

        assert expected_keepD == keepD
        assert expected_trimD == trimD
        assert logArr == []

class TestJoinKeepDAndTrimD(object):

    def test_join_keepD_and_trimD(self):
        ## set up
        keepD = {
            '1': ['A', '-', 'G', 'T', 'A', 'T'], 
            '2': ['A', '-', 'G', '-', 'A', 'T'], 
            '3': ['A', '-', 'G', '-', 'T', 'A'], 
            '4': ['A', 'G', 'A', '-', 'T', 'A'], 
            '5': ['A', 'C', 'a', '-', 'T', '-']
            }
        trimD = {'1': [], '2': [], '3': [], '4': [], '5': []}
        logArr = [
            ['1', 'keep', 'nPI', '0.0'],
            ['2', 'keep', 'nPI', '0.6'],
            ['3', 'keep', 'PI', '0.0'],
            ['4', 'keep', 'nPI', '0.8'],
            ['5', 'keep', 'PI', '0.0'],
            ['6', 'keep', 'PI', '0.2']
            ]

        ## execution
        keepD, trimD = join_keepD_and_trimD(keepD, trimD)

        ## check results
        expected_keepD = {'1': 'A-GTAT',
            '2': 'A-G-AT',
            '3': 'A-G-TA',
            '4': 'AGA-TA',
            '5': 'ACa-T-'
            }
        expected_trimD = {'1': '', '2': '', '3': '', '4': '', '5': ''}

        assert expected_keepD == keepD
        assert expected_trimD == trimD

## TODO: finish writing
class TestWriteTrimD(object):

    def test_write_trimD(self, mocker):
        ## set up
        trimD = {
            '1':['A'],
            '2':['A'],
            '3':['A'],
            '4':['A'],
            '5':['A']
            }
        outFile = 'output_file_name.fa'
        outFileFormat = 'fasta'
        fake_msa = MultipleSeqAlignment(
            [SeqRecord(seq=Seq("['A']"),
                id='1',
                name='<unknown name>',
                description='',
                dbxrefs=[]), 
            SeqRecord(seq=Seq("['A']"),
                id='2',
                name='<unknown name>',
                description='',
                dbxrefs=[]),
                SeqRecord(seq=Seq("['A']"),
                id='3',
                name='<unknown name>',
                description='', 
                dbxrefs=[]),
            SeqRecord(seq=Seq("['A']"),
                id='4',
                name='<unknown name>',
                description='',
                dbxrefs=[]),
            SeqRecord(seq=Seq("['A']"),
                id='5',
                name='<unknown name>',
                description='',
                dbxrefs=[])]
                )
        fake_msa.__repr__ = 'a3c184c'
        mock_msa = mocker.patch('Bio.Align.MultipleSeqAlignment')
        mock_msa.return_value = (f"{here.parent}/examples/simple.fa", "fasta")
        mock_write = mocker.patch("Bio.SeqIO.write")

        ## execution
        write_trimD(trimD, outFileFormat, outFile)

        ## check results
        mock_write.assert_called_once_with(fake_msa, f"${outFile}.complement", outFileFormat)

        
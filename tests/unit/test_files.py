import pytest
from Bio import AlignIO
from pathlib import Path

here = Path(__file__)

from clipkit.files import automatic_file_type_determination, FileFormat

class TestAutomaticFileTypeDetermination(object):

    def test_automatic_file_type_determination(self):
        # setup
        in_file = f"{here.parent}/examples/simple.fa"

        # execution
        alignment, in_file_format = automatic_file_type_determination(in_file)

        # check results
        assert in_file_format == FileFormat.fasta
        assert alignment.get_alignment_length() == 6

    def test_automatic_file_type_determination_raises_error_when_file_not_known(self, mocker):
        in_file = ''
        mocker.patch('clipkit.files.AlignIO.read', side_effect=ValueError())

        with pytest.raises(Exception):
            automatic_file_type_determination(in_file)

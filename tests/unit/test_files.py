import pytest
from Bio import AlignIO
from pathlib import Path

here = Path(__file__)

from clipkit.files import get_alignment_and_format, FileFormat

class TestAutomaticFileTypeDetermination(object):

    def test_get_alignment_and_format_when_format_is_provided(self):
        # setup
        in_file = f"{here.parent}/examples/simple.fa"
        file_format = FileFormat.fasta

        # execution
        alignment, in_file_format = get_alignment_and_format(in_file, file_format)

        # check results
        assert in_file_format == FileFormat.fasta
        assert alignment.get_alignment_length() == 6

    def test_get_alignment_and_format_when_format_is_not_provided(self):
        # setup
        in_file = f"{here.parent}/examples/simple.fa"
        file_format = None

        # execution
        alignment, in_file_format = get_alignment_and_format(in_file, file_format)

        # check results
        assert in_file_format == FileFormat.fasta
        assert alignment.get_alignment_length() == 6

    def test_get_alignment_and_format_raises_error_when_file_not_known(self, mocker):
        in_file = ''
        file_format = None
        mocker.patch('clipkit.files.AlignIO.read', side_effect=ValueError())

        with pytest.raises(Exception):
            get_alignment_and_format(in_file, file_format)

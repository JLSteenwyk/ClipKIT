import pytest
from pathlib import Path

from Bio import AlignIO
from clipkit.files import get_alignment_and_format, FileFormat

here = Path(__file__)


class TestAutomaticFileTypeDetermination(object):
    def test_get_alignment_and_format_when_format_is_provided(self):
        # setup
        in_file = f"{here.parent}/examples/simple.fa"
        file_format = "fasta"

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
        in_file = ""
        file_format = None
        mocker.patch("clipkit.files.AlignIO.read", side_effect=ValueError())

        with pytest.raises(Exception) as excinfo:
            get_alignment_and_format(in_file, file_format)
        assert "No such file or directory" in str(excinfo.value)

    def test_get_alignment_and_format_raises_error_when_detection_fails(self, mocker):
        in_file = f"{here.parent}/examples/simple.fa"
        file_format = None
        mocker.patch("clipkit.files.AlignIO.read", side_effect=ValueError())

        with pytest.raises(Exception) as excinfo:
            get_alignment_and_format(in_file, file_format)
        assert "File could not be read" in str(excinfo.value)

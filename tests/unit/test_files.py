import pytest
from Bio import AlignIO
from pathlib import Path

here = Path(__file__)

from clipkit.files import automatic_file_type_determination

class TestAutomaticFileTypeDetermination(object):

    def test_automatic_file_type_determination(self):
        # setup
        in_file = f"{here.parent}/examples/simple.fa"

        # execution
        alignment, in_file_format = automatic_file_type_determination(in_file)

        # check results
        assert in_file_format == "fasta"
        assert alignment.get_alignment_length() == 6
        
# ASK THOMAS HOW TO TEST A PRINT MESSAGE


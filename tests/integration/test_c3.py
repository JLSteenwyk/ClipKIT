import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode
from clipkit.settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS

here = Path(__file__)


@pytest.mark.integration
class TestC3Out(object):
    def test_simple_c3(self):
        """
        test codon
        usage: clipkit simple.fa c3
        """
        output_file = "output/simple.fa_c3"

        kwargs = dict(
            input_file=f"{here.parent}/samples/simple.fa",
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=None,
            mode=TrimmingMode.c3,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )

        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.fa_c3", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

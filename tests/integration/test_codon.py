import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode
from clipkit.settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS

here = Path(__file__)


@pytest.mark.integration
class TestCodonOut(object):
    def test_simple_codon(self):
        """
        test codon
        usage: clipkit simple.fa -co
        """
        output_file = "output/simple.fa_gappy_codon"

        kwargs = dict(
            input_file=f"{here.parent}/samples/simple.fa",
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            codon=True,
            gaps=0.8,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
            ends_only=False,
        )

        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.fa_gappy_codon", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_simple_codon_all_trimmed(self):
        """
        test codon
        usage: clipkit simple.fa -co
        """
        output_file = "output/simple.fa_gappy_codon"
        kwargs = dict(
            input_file=f"{here.parent}/samples/simple.fa",
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            codon=True,
            gaps=0.1,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
            ends_only=False,
        )

        execute(**kwargs)

        expected_content = ">1\n>2\n>3\n>4\n>5\n"

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

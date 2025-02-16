import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode
from clipkit.settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS

here = Path(__file__)


@pytest.mark.integration
class TestEndsOnly(object):
    def test_ends_only(self):
        """
        test smart-gap with simple_extra_gap_trim_ends_example.fa
        usage: clipkit simple_extra_gap_trim_ends_example.fa
        """
        input_file = f"{here.parent}/samples/simple_extra_gap_trim_ends_example.fa"
        output_file = "output/simple_extra_gap_trim_ends_example.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.8,
            mode=TrimmingMode.smart_gap,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
            ends_only=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple_extra_gap_trim_ends_example.fa.clipkit", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

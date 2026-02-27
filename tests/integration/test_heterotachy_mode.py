import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.modes import TrimmingMode
from clipkit.settings import DEFAULT_NT_GAP_CHARS

here = Path(__file__)


@pytest.mark.integration
class TestHeterotachyMode(object):
    def test_simple_alignment(self):
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.fa_heterotachy"
        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.8,
            mode=TrimmingMode.heterotachy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
            ends_only=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.fa_heterotachy", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

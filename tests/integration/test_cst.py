import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode
from clipkit.settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS

here = Path(__file__)


@pytest.mark.integration
class TestCSTOut(object):
    @pytest.mark.parametrize(
        "aux_file, expected_result_file",
        [
            ("samples/cst_ex0.txt", "expected/simple_cst_ex0.fa"),
            ("samples/cst_ex1.txt", "expected/simple_cst_ex1.fa"),
            ("samples/cst_ex2.txt", "expected/simple_cst_ex2.fa"),
            ("samples/cst_ex3.txt", "expected/simple_cst_ex3.fa"),
        ],
    )
    def test_cst(self, aux_file, expected_result_file):
        """
        test custom site trimming
        usage: clipkit simple.fa -m cst -a tests/integration/samples/cst_ex0.txt
        """
        output_file = "output/simple.fa_cst"

        kwargs = dict(
            input_file=f"{here.parent}/samples/simple.fa",
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=None,
            mode=TrimmingMode.cst,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
            auxiliary_file=f"{here.parent}/{aux_file}",
        )

        execute(**kwargs)

        with open(f"{here.parent}/{expected_result_file}", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

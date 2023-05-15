import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)


@pytest.mark.integration
class TestOutLog(object):
    def test_simple(self):
        """
        test output in clustal format
        usage: clipkit simple.fa -l
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=True,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.fa.clipkit.log", "r") as expected:
            expected_content = expected.read()

        with open(f"{output_file}.log", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_simple_long_description(self):
        """
        test output in clustal format
        usage: clipkit simple_long_description.fa -l
        """
        input_file = f"{here.parent}/samples/simple_long_description.fa"
        output_file = "output/simple_long_description.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple_long_description.fa.clipkit.log", "r"
        ) as expected:
            expected_content = expected.read()

        with open(f"{output_file}.log", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_codon_aln(self):
        """
        test output in clustal format
        usage: clipkit 12_YIL115C_Anc_2.253_codon_aln.fasta -l
        """
        input_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_codon_aln.fasta"
        output_file = "output/12_YIL115C_Anc_2.253_codon_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_codon_aln.fasta.clipkit.log",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(f"{output_file}.log", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG091N44M8_nt(self):
        """
        test output in clustal format
        usage: clipkit EOG091N44M8_nt.fa -l
        """
        input_file = f"{here.parent}/samples/EOG091N44M8_nt.fa"
        output_file = "output/EOG091N44M8_nt.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="fasta",
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG091N44M8_nt.fa.clipkit.log", "r"
        ) as expected:
            expected_content = expected.read()

        with open(f"{output_file}.log", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

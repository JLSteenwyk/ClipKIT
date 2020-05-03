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
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.fa.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=True,
        )

        with open(f"{here.parent}/expected/simple.fa.clipkit.log", "r") as expected:
            expected_content = expected.read()

        with open(f"{out_file}.log", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_codon_aln(self):
        """
        test output in clustal format
        usage: clipkit 12_YIL115C_Anc_2.253_codon_aln.fasta -l
        """
        in_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_codon_aln.fasta"
        out_file = "output/12_YIL115C_Anc_2.253_codon_aln.fasta.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=True,
        )

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_codon_aln.fasta.clipkit.log",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(f"{out_file}.log", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG091N44M8_nt(self):
        """
        test output in clustal format
        usage: clipkit EOG091N44M8_nt.fa -l
        """
        in_file = f"{here.parent}/samples/EOG091N44M8_nt.fa"
        out_file = "output/EOG091N44M8_nt.fa.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=True,
        )

        with open(
            f"{here.parent}/expected/EOG091N44M8_nt.fa.clipkit.log", "r"
        ) as expected:
            expected_content = expected.read()

        with open(f"{out_file}.log", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

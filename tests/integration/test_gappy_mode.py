import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)

@pytest.mark.integration
class TestGappyMode(object):
    def test_simple_no_change(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simpla.fa.clipkit"
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
        )

        with open(in_file, 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2(self):
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
        )

        with open(f"{here.parent}/expected/12_YIL115C_Anc_2.253_codon_aln.fasta", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

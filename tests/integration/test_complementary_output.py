import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)

@pytest.mark.integration
class TestComplementOut(object):
    # test complementary output file with a simple case
    # usage: clipkit simple.fa -c
    def test_simple_complement(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.fa_gappy"
        complement_out_file = f"{out_file}.complement"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=True,
            mode=TrimmingMode.gappy,
        )

        with open(f"{here.parent}/expected/simple.fa_gappy.complement", 'r') as expected:
            expected_content = expected.read()

        with open(complement_out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test complementary output file for amino acid yeast sequences
    # usage: clipkit 12_YIL115C_Anc_2.253_aa_aln.fasta -c
    def test_12_YIL115C_Anc_2_253_aa_aln_complement(self):
        in_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_aa_aln.fasta"
        out_file = "output/12_YIL115C_Anc_2.253_aa_aln.fasta_gappy"
        complement_out_file = f"{out_file}.complement"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=True,
            mode=TrimmingMode.gappy,
        )

        with open(f"{here.parent}/expected/12_YIL115C_Anc_2.253_aa_aln.fasta_gappy.complement", 'r') as expected:
            expected_content = expected.read()

        with open(complement_out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test complementary output file for amino acid Penicillium sequences
    # usage: clipkit EOG091N44M8_aa.fa -c
    def test_EOG091N44M8_aa_complement(self):
        in_file = f"{here.parent}/samples/EOG091N44M8_aa.fa"
        out_file = "output/EOG091N44M8_aa.fa_gappy"
        complement_out_file = f"{out_file}.complement"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=True,
            mode=TrimmingMode.gappy,
        )

        with open(f"{here.parent}/expected/EOG091N44M8_aa.fa_gappy.complement", 'r') as expected:
            expected_content = expected.read()

        with open(complement_out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content
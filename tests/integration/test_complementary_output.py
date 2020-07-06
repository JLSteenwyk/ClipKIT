import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)


@pytest.mark.integration
class TestComplementOut(object):
    def test_simple_complement(self):
        """
        test complementary output file with a simple case
        usage: clipkit simple.fa -c
        """
        output_file = "output/simple.fa_gappy"
        complement_out_file = f"{output_file}.complement"

        kwargs = dict(
            input_file=f"{here.parent}/samples/simple.fa",
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            complement=True,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )

        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple.fa_gappy.complement", "r"
        ) as expected:
            expected_content = expected.read()

        with open(complement_out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_aa_aln_complement(self):
        """
        test complementary output file for amino acid yeast sequences
        usage: clipkit 12_YIL115C_Anc_2.253_aa_aln.fasta -c
        """
        output_file = "output/12_YIL115C_Anc_2.253_aa_aln.fasta_gappy"
        complement_out_file = f"{output_file}.complement"

        kwargs = dict(
            input_file=f"{here.parent}/samples/12_YIL115C_Anc_2.253_aa_aln.fasta",
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            complement=True,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_aa_aln.fasta_gappy.complement",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(complement_out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG091N44M8_aa_complement(self):
        """
        test complementary output file for amino acid Penicillium sequences
        usage: clipkit EOG091N44M8_aa.fa -c
        """
        output_file = "output/EOG091N44M8_aa.fa_gappy"
        complement_out_file = f"{output_file}.complement"

        kwargs = dict(
            input_file=f"{here.parent}/samples/EOG091N44M8_aa.fa",
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            complement=True,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG091N44M8_aa.fa_gappy.complement", "r"
        ) as expected:
            expected_content = expected.read()

        with open(complement_out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

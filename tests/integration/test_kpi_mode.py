import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)


@pytest.mark.integration
class TestKpiMode(object):
    def test_simple(self):
        """
        usage: clipkit simple.fa -kpi
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simpla.fa.TestKpiMode_test_simple.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.fa_kpi", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_simple_long_description(self):
        """
        usage: clipkit simple.fa -kpi
        """
        input_file = f"{here.parent}/samples/simple_long_description.fa"
        output_file = "output/simple_long_description.fa.TestKpiMode_test_simple.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple_long_description.fa_kpi", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_codon_aln(self):
        """
        test kpi with codon alignment of yeast sequences
        usage: clipkit 12_YIL115C_Anc_2.253_codon_aln.fasta -m kpi
        """
        input_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_codon_aln.fasta"
        output_file = "output/12_YIL115C_Anc_2.253_codon_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_codon_aln.fasta_kpi", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_aa_aln(self):
        """
        test kpi with amino acid alignment of yeast sequences
        usage: clipkit 12_YIL115C_Anc_2.253_aa_aln.fasta -m kpi
        """
        input_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_aa_aln.fasta"
        output_file = "output/12_YIL115C_Anc_2.253_aa_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_aa_aln.fasta_kpi", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_24_ENSG00000163519_aa_aln(self):
        """
        test kpi with amino acid alignment of mammalian sequences
        usage: clipkit 24_ENSG00000163519_aa_aln.fasta -m kpi
        """
        input_file = f"{here.parent}/samples/24_ENSG00000163519_aa_aln.fasta"
        output_file = "output/24_ENSG00000163519_aa_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/24_ENSG00000163519_aa_aln.fasta_kpi", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_24_ENSG00000163519_codon_aln(self):
        """
        test kpi with codon alignment of mammalian sequences
        usage: clipkit 24_ENSG00000163519_codon_aln.fasta -m kpi
        """
        input_file = f"{here.parent}/samples/24_ENSG00000163519_codon_aln.fasta"
        output_file = "output/24_ENSG00000163519_codon_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/24_ENSG00000163519_codon_aln.fasta_kpi", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG091N44M8_aa(self):
        """
        test kpi with amino acid alignment of Penicillium sequences
        usage: clipkit EOG091N44M8_aa.fa -m kpi
        """
        input_file = f"{here.parent}/samples/EOG091N44M8_aa.fa"
        output_file = "output/EOG091N44M8_aa.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/EOG091N44M8_aa.fa_kpi", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG091N44M8_nt(self):
        """
        test kpi with nucleotide alignment of Penicillium sequences
        usage: clipkit EOG091N44M8_nt.fa -m kpi
        """
        input_file = f"{here.parent}/samples/EOG091N44M8_nt.fa"
        output_file = "output/EOG091N44M8_nt.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format='fasta',
            output_file_format='fasta',
            sequence_type=None,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/EOG091N44M8_nt.fa_kpi", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # @pytest.mark.slow
    # def test_EOG092C4VOX_aa(self):
    #     """
    #     test kpi with amino alignment of fungal sequences
    #     usage: clipkit EOG092C4VOX_aa_aln.fasta -m kpi
    #     """
    #     input_file = f"{here.parent}/samples/EOG092C4VOX_aa_aln.fasta"
    #     output_file = "output/EOG092C4VOX_aa_aln.fasta.clipkit"

    #     kwargs = dict(
    #         input_file=input_file,
    #         output_file=output_file,
    #         input_file_format='fasta',
    #         output_file_format='fasta',
    #         complement=False,
    #         gaps=0.9,
    #         mode=TrimmingMode.kpi,
    #         use_log=False,
    #     )
    #     execute(**kwargs)

    #     with open(
    #         f"{here.parent}/expected/EOG092C4VOX_aa_aln.fasta_kpi", "r"
    #     ) as expected:
    #         expected_content = expected.read()

    #     with open(output_file, "r") as out_file:
    #         output_content = out_file.read()

    #     assert expected_content == output_content

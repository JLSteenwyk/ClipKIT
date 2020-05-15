import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)


@pytest.mark.integration
class TestKpiGappyMode(object):
    def test_simple(self):
        """
        usage: clipkit simple.fa -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simpla.fa.TestKpiGappyMode_test_simple.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.fa_kpi_gappy", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_codon_aln(self):
        """
        test kpi_gappy with codon alignment of yeast sequences
        usage: clipkit 12_YIL115C_Anc_2.253_codon_aln.fasta -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_codon_aln.fasta"
        output_file = "output/12_YIL115C_Anc_2.253_codon_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_codon_aln.fasta_kpi_gappy",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_aa_aln(self):
        """
        test kpi_gappy with amino acid alignment of yeast sequences
        usage: clipkit 12_YIL115C_Anc_2.253_aa_aln.fasta -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_aa_aln.fasta"
        output_file = "output/12_YIL115C_Anc_2.253_aa_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_aa_aln.fasta_kpi_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_24_ENSG00000163519_aa_aln(self):
        """
        test kpi_gappy with amino acid alignment of mammalian sequences
        usage: clipkit 24_ENSG00000163519_aa_aln.fasta -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/24_ENSG00000163519_aa_aln.fasta"
        output_file = "output/24_ENSG00000163519_aa_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/24_ENSG00000163519_aa_aln.fasta_kpi_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_24_ENSG00000163519_codon_aln(self):
        """
        test kpi_gappy with codon alignment of mammalian sequences
        usage: clipkit 24_ENSG00000163519_codon_aln.fasta -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/24_ENSG00000163519_codon_aln.fasta"
        output_file = "output/24_ENSG00000163519_codon_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/24_ENSG00000163519_codon_aln.fasta_kpi_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG091N44M8_aa(self):
        """
        test kpi_gappy with amino acid alignment of Penicillium sequences
        usage: clipkit EOG091N44M8_aa.fa -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/EOG091N44M8_aa.fa"
        output_file = "output/EOG091N44M8_aa.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG091N44M8_aa.fa_kpi_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG091N44M8_nt(self):
        """
        test kpi_gappy with nucleotide alignment of Penicillium sequences
        usage: clipkit EOG091N44M8_nt.fa -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/EOG091N44M8_nt.fa"
        output_file = "output/EOG091N44M8_nt.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG091N44M8_nt.fa_kpi_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    @pytest.mark.slow
    def test_EOG092C0CZK_aa(self):
        """
        test kpi_gappy with amino alignment of fungal sequences
        usage: clipkit EOG092C0CZK_aa_aln.fasta -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/EOG092C0CZK_aa_aln.fasta"
        output_file = "output/EOG092C0CZK_aa_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG092C0CZK_aa_aln.fasta_kpi_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG092C4VOX_aa(self):
        """
        test gappy with amino alignment of fungal sequences
        usage: clipkit EOG092C4VOX_aa_aln.fasta -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/EOG092C4VOX_aa_aln.fasta"
        output_file = "output/EOG092C4VOX_aa_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG092C4VOX_aa_aln.fasta_kpi_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content


@pytest.mark.integration
class TestKPIGappyModeCustomGapsParameter(object):
    def test_simple(self):
        """
        test kpi_gappy with a custom gaps parameter
        usage: clipkit simple.fa -g 0.2 -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simpla.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.2,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple.fa_kpi_gappy_gaps_set_to_0.2", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_codon_aln(self):
        """
        test kpi_gappy with codon alignment of yeast sequences
        usage: clipkit 12_YIL115C_Anc_2.253_codon_aln.fasta -g 0.3 -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_codon_aln.fasta"
        output_file = "output/12_YIL115C_Anc_2.253_codon_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.3,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_codon_aln.fasta_kpi_gappy_custom_gaps",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_24_ENSG00000163519_codon_aln(self):
        """
        test kpi_gappy with codon alignment of mammalian sequences
        usage: clipkit 24_ENSG00000163519_codon_aln.fasta -g .4 -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/24_ENSG00000163519_codon_aln.fasta"
        output_file = "output/24_ENSG00000163519_codon_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.4,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/24_ENSG00000163519_codon_aln.fasta_kpi_gappy_custom_gaps",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG091N44M8_nt(self):
        """
        test kpi_gappy with nucleotide alignment of Penicillium sequences
        usage: clipkit EOG091N44M8_nt.fa -g .1 -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/EOG091N44M8_nt.fa"
        output_file = "output/EOG091N44M8_nt.fa.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.1,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG091N44M8_nt.fa_kpi_gappy_custom_gaps", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    @pytest.mark.slow
    def test_EOG092C0CZK_aa(self):
        """
        test kpi_gappy with amino alignment of fungal sequences
        usage: clipkit EOG092C0CZK_aa_aln.fasta -g .5 -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/EOG092C0CZK_aa_aln.fasta"
        output_file = "output/EOG092C0CZK_aa_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.5,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG092C0CZK_aa_aln.fasta_kpi_gappy_custom_gaps",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_EOG092C4VOX_aa(self):
        """
        test kpi_gappy with amino alignment of fungal sequences
        usage: clipkit EOG092C4VOX_aa_aln.fasta -g .25 -m kpi-gappy
        """
        input_file = f"{here.parent}/samples/EOG092C4VOX_aa_aln.fasta"
        output_file = "output/EOG092C4VOX_aa_aln.fasta.clipkit"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.fasta,
            complement=False,
            gaps=0.25,
            mode=TrimmingMode.kpi_gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/EOG092C4VOX_aa_aln.fasta_kpi_gappy_custom_gaps",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

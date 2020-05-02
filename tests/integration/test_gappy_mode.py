import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)


@pytest.mark.integration
class TestGappyMode(object):
    def test_simple_no_change(self):
        """
        test gappy where no changes are expected in the resulting 
        output alignment.
        usage: clipkit simple.fa
        """
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
            use_log=False,
        )

        with open(in_file, "r") as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_codon_aln(self):
        """
        test gappy with codon alignment of yeast sequences
        usage: clipkit 12_YIL115C_Anc_2.253_codon_aln.fasta
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
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_codon_aln.fasta_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_12_YIL115C_Anc_2_253_aa_aln(self):
        """
        test gappy with amino acid alignment of yeast sequences
        usage: clipkit 12_YIL115C_Anc_2.253_aa_aln.fasta 
        """
        in_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_aa_aln.fasta"
        out_file = "output/12_YIL115C_Anc_2.253_aa_aln.fasta.clipkit"
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
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_aa_aln.fasta_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_24_ENSG00000163519_aa_aln(self):
        """
        test gappy with amino acid alignment of mammalian sequences
        usage: clipkit 24_ENSG00000163519_aa_aln.fasta 
        """
        in_file = f"{here.parent}/samples/24_ENSG00000163519_aa_aln.fasta"
        out_file = "output/24_ENSG00000163519_aa_aln.fasta.clipkit"
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
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/24_ENSG00000163519_aa_aln.fasta_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with codon alignment of mammalian sequences
    # usage: clipkit 24_ENSG00000163519_codon_aln.fasta
    def test_24_ENSG00000163519_codon_aln(self):
        in_file = f"{here.parent}/samples/24_ENSG00000163519_codon_aln.fasta"
        out_file = "output/24_ENSG00000163519_codon_aln.fasta.clipkit"
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
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/24_ENSG00000163519_codon_aln.fasta_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with amino acid alignment of Penicillium sequences
    # usage: clipkit EOG091N44M8_aa.fa
    def test_EOG091N44M8_aa(self):
        in_file = f"{here.parent}/samples/EOG091N44M8_aa.fa"
        out_file = "output/EOG091N44M8_aa.fa.clipkit"
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
            use_log=False,
        )

        with open(f"{here.parent}/expected/EOG091N44M8_aa.fa_gappy", "r") as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with nucleotide alignment of Penicillium sequences
    # usage: clipkit EOG091N44M8_nt.fa
    def test_EOG091N44M8_nt(self):
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
            use_log=False,
        )

        with open(f"{here.parent}/expected/EOG091N44M8_nt.fa_gappy", "r") as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with amino alignment of fungal sequences
    # usage: clipkit EOG092C0CZK_aa_aln.fasta
    def test_EOG092C0CZK_aa(self):
        in_file = f"{here.parent}/samples/EOG092C0CZK_aa_aln.fasta"
        out_file = "output/EOG092C0CZK_aa_aln.fasta.clipkit"
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
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/EOG092C0CZK_aa_aln.fasta_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with amino alignment of fungal sequences
    # usage: clipkit EOG092C4VOX_aa_aln.fasta
    def test_EOG092C4VOX_aa(self):
        in_file = f"{here.parent}/samples/EOG092C4VOX_aa_aln.fasta"
        out_file = "output/EOG092C4VOX_aa_aln.fasta.clipkit"
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
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/EOG092C4VOX_aa_aln.fasta_gappy", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content


@pytest.mark.integration
class TestGappyModeCustomGapsParameter(object):
    # test gappy with a custom gaps parameter
    # usage: clipkit simple.fa -g 0.2
    def test_simple(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simpla.fa.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.2,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/simple.fa_gappy_gaps_set_to_0.2", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with codon alignment of yeast sequences
    # usage: clipkit 12_YIL115C_Anc_2.253_codon_aln.fasta -g 0.3
    def test_12_YIL115C_Anc_2_253_codon_aln(self):
        in_file = f"{here.parent}/samples/12_YIL115C_Anc_2.253_codon_aln.fasta"
        out_file = "output/12_YIL115C_Anc_2.253_codon_aln.fasta.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.3,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/12_YIL115C_Anc_2.253_codon_aln.fasta_gappy_custom_gaps",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with codon alignment of mammalian sequences
    # usage: clipkit 24_ENSG00000163519_codon_aln.fasta -g .4
    def test_24_ENSG00000163519_codon_aln(self):
        in_file = f"{here.parent}/samples/24_ENSG00000163519_codon_aln.fasta"
        out_file = "output/24_ENSG00000163519_codon_aln.fasta.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.4,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/24_ENSG00000163519_codon_aln.fasta_gappy_custom_gaps",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with nucleotide alignment of Penicillium sequences
    # usage: clipkit EOG091N44M8_nt.fa -g .1
    def test_EOG091N44M8_nt(self):
        in_file = f"{here.parent}/samples/EOG091N44M8_nt.fa"
        out_file = "output/EOG091N44M8_nt.fa.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.1,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/EOG091N44M8_nt.fa_gappy_custom_gaps", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with amino alignment of fungal sequences
    # usage: clipkit EOG092C0CZK_aa_aln.fasta -g .5
    def test_EOG092C0CZK_aa(self):
        in_file = f"{here.parent}/samples/EOG092C0CZK_aa_aln.fasta"
        out_file = "output/EOG092C0CZK_aa_aln.fasta.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.5,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/EOG092C0CZK_aa_aln.fasta_gappy_custom_gaps", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test gappy with amino alignment of fungal sequences
    # usage: clipkit EOG092C4VOX_aa_aln.fasta -g .25
    def test_EOG092C4VOX_aa(self):
        in_file = f"{here.parent}/samples/EOG092C4VOX_aa_aln.fasta"
        out_file = "output/EOG092C4VOX_aa_aln.fasta.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.25,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False,
        )

        with open(
            f"{here.parent}/expected/EOG092C4VOX_aa_aln.fasta_gappy_custom_gaps", "r"
        ) as expected:
            expected_content = expected.read()

        with open(out_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

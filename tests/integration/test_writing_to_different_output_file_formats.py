import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode
from clipkit.settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS

here = Path(__file__)


@pytest.mark.integration
class TestOutFormats(object):
    def test_clustal(self):
        """
        test output in clustal format
        usage: clipkit simple.fa -of clustal
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.clustal"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="clustal",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.clustal", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_maf(self):
        """
        test output in maf format
        usage: clipkit simple.fa -of maf
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.maf"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="maf",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.maf", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_mauve(self):
        """
        test output in mauve format
        usage: clipkit simple.fa -of mauve
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.mauve"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="mauve",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.mauve", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_phylip(self):
        """
        test output in phylip format
        usage: clipkit simple.fa -of phylip
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.phylip"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="phylip",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.phylip", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_phylip_sequential(self):
        """
        test output in phylip-sequential format
        usage: clipkit simple.fa -of phylip-sequential
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.phylip-sequential"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="phylip_sequential",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.phylip-sequential", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_phylip_relaxed(self):
        """
        test output in phylip-relaxed format
        usage: clipkit simple.fa -of phylip-relaxed
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.phylip-relaxed"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="phylip_relaxed",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.phylip-relaxed", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_stockholm(self):
        """
        test output in stockholm format
        usage: clipkit simple.fa -of stockholm
        """
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.stockholm"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="stockholm",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.stockholm", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_clustal_long_description(self):
        """
        test output in clustal format
        usage: clipkit simple_long_description.fa -of clustal
        """
        input_file = f"{here.parent}/samples/simple_long_description.fa"
        output_file = "output/simple_long_description.clustal"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="clustal",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple_long_description.clustal", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_maf_long_description(self):
        """
        test output in maf format
        usage: clipkit simple_long_description.fa -of maf
        """
        input_file = f"{here.parent}/samples/simple_long_description.fa"
        output_file = "output/simple_long_description.maf"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="maf",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple_long_description.maf", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_mauve_long_description(self):
        """
        test output in mauve format
        usage: clipkit simple_long_description.fa -of mauve
        """
        input_file = f"{here.parent}/samples/simple_long_description.fa"
        output_file = "output/simple_long_description.mauve"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="mauve",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple_long_description.mauve", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_phylip_long_description(self):
        """
        test output in phylip format
        usage: clipkit simple_long_description.fa -of phylip
        """
        input_file = f"{here.parent}/samples/simple_long_description.fa"
        output_file = "output/simple_long_description.phylip"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="phylip",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple_long_description.phylip", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_phylip_sequential_long_description(self):
        """
        test output in phylip-sequential format
        usage: clipkit simple_long_description.fa -of phylip-sequential
        """
        input_file = f"{here.parent}/samples/simple_long_description.fa"
        output_file = "output/simple_long_description.phylip-sequential"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="phylip_sequential",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple_long_description.phylip-sequential", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_stockholm_long_description(self):
        """
        test output in stockholm format
        usage: clipkit simple_long_description.fa -of stockholm
        """
        input_file = f"{here.parent}/samples/simple_long_description.fa"
        output_file = "output/simple_long_description.stockholm"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format="fasta",
            output_file_format="stockholm",
            sequence_type=None,
            complement=False,
            codon=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
            gap_characters=DEFAULT_NT_GAP_CHARS,
            quiet=True,
        )
        execute(**kwargs)

        with open(
            f"{here.parent}/expected/simple_long_description.stockholm", "r"
        ) as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

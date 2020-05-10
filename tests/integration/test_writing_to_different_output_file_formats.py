import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)


@pytest.mark.integration
class TestOutFormats(object):
    # test output in clustal format
    # usage: clipkit simple.fa -of clustal
    def test_clustal(self):
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.clustal"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.clustal,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.clustal", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in maf format
    # usage: clipkit simple.fa -of maf
    def test_maf(self):
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.maf"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.maf,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.maf", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in mauve format
    # usage: clipkit simple.fa -of mauve
    def test_mauve(self):
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.mauve"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.mauve,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.mauve", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in phylip format
    # usage: clipkit simple.fa -of phylip
    def test_phylip(self):
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.phylip"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.phylip,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.phylip", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in phylip-sequential format
    # usage: clipkit simple.fa -of phylip-sequential
    def test_phylip_sequential(self):
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.phylip-sequential"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.phylip_seq,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.phylip-sequential", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in phylip-relaxed format
    # usage: clipkit simple.fa -of phylip-relaxed
    def test_phylip_relaxed(self):
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.phylip-relaxed"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.phylip_rel

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.phylip_rel,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.phylip-relaxed", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in stockholm format
    # usage: clipkit simple.fa -of stockholm
    def test_stockholm(self):
        input_file = f"{here.parent}/samples/simple.fa"
        output_file = "output/simple.stockholm"

        kwargs = dict(
            input_file=input_file,
            output_file=output_file,
            input_file_format=FileFormat.fasta,
            output_file_format=FileFormat.stockholm,
            complement=False,
            gaps=0.9,
            mode=TrimmingMode.gappy,
            use_log=False,
        )
        execute(**kwargs)

        with open(f"{here.parent}/expected/simple.stockholm", "r") as expected:
            expected_content = expected.read()

        with open(output_file, "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

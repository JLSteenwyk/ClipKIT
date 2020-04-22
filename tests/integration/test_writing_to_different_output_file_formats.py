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
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.clustal"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.clustal

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False
        )

        with open(f"{here.parent}/expected/simple.clustal", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in maf format
    # usage: clipkit simple.fa -of maf
    def test_maf(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.maf"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.maf

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False
        )

        with open(f"{here.parent}/expected/simple.maf", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in mauve format
    # usage: clipkit simple.fa -of mauve
    def test_mauve(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.mauve"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.mauve

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False
        )

        with open(f"{here.parent}/expected/simple.mauve", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in phylip format
    # usage: clipkit simple.fa -of phylip
    def test_phylip(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.phylip"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.phylip

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False
        )

        with open(f"{here.parent}/expected/simple.phylip", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in phylip-sequential format
    # usage: clipkit simple.fa -of phylip-sequential
    def test_phylip_sequential(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.phylip-sequential"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.phylip_seq

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False
        )

        with open(f"{here.parent}/expected/simple.phylip-sequential", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    # test output in phylip-relaxed format
    # usage: clipkit simple.fa -of phylip-relaxed
    def test_phylip_relaxed(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.phylip-relaxed"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.phylip_rel

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False
        )

        with open(f"{here.parent}/expected/simple.phylip-relaxed", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content
    
    # test output in stockholm format
    # usage: clipkit simple.fa -of stockholm
    def test_stockholm(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simple.stockholm"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.stockholm

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.gappy,
            use_log=False
        )

        with open(f"{here.parent}/expected/simple.stockholm", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content
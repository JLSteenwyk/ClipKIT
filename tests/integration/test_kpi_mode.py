import pytest
from pathlib import Path

from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode

here = Path(__file__)

@pytest.mark.integration
class TestKpiMode(object):
    def test_simple(self):
        in_file = f"{here.parent}/samples/simple.fa"
        out_file = "output/simpla.fa.TestKpiMode_test_simple.clipkit"
        in_file_format = FileFormat.fasta
        out_file_format = FileFormat.fasta

        execute(
            in_file,
            out_file,
            in_file_format,
            out_file_format,
            gaps=0.9,
            complement=False,
            mode=TrimmingMode.kpi,
        )

        with open(f"{here.parent}/expected/simple.fa_kpi", 'r') as expected:
            expected_content = expected.read()

        with open(out_file, 'r') as out_file:
            output_content = out_file.read()

        assert expected_content == output_content



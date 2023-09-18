import pytest

from Bio.Align import MultipleSeqAlignment
from clipkit import clipkit
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode
from clipkit.msa import MSA

@pytest.mark.integration
class TestApiInvocation(object):
    def test_input_file(self):
        trim_run, stats = clipkit(
            input_file_path="tests/integration/samples/simple.fa",
            mode=TrimmingMode.gappy,
            gaps=0.3,
            sequence_type="nt",
        )

        assert stats.summary == {
            "alignment_length": 6,
            "output_length": 4,
            "trimmed_length": 2,
            "trimmed_percentage": 33.333
        }
        assert isinstance(trim_run.version, str)
        assert isinstance(trim_run.trimmed, MultipleSeqAlignment)

    def test_raw_alignment(self):
        trim_run, stats = clipkit(
            raw_alignment=">1\nA-GTAT\n>2\nA-G-AT\n>3\nA-G-TA\n>4\nAGA-TA\n>5\nACa-T-\n",
            mode=TrimmingMode.smart_gap,
            gaps=None,
            sequence_type="nt",
        )
        assert stats.summary == {
            "alignment_length": 6,
            "output_length": 5,
            "trimmed_length": 1,
            "trimmed_percentage": 16.667
        }
        assert isinstance(trim_run.version, str)

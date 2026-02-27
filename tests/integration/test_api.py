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
            "trimmed_percentage": 33.333,
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
            "trimmed_percentage": 16.667,
        }
        assert isinstance(trim_run.version, str)
        
    def test_codon_setting(self):
        trim_run, stats = clipkit(
            raw_alignment=">1\nA-GTAT\n>2\nA-G-AT\n>3\nA-G-TA\n>4\nAGA-TA\n>5\nACa-T-\n",
            mode=TrimmingMode.smart_gap,
            gaps=None,
            codon=True,
            sequence_type="nt",
        )
        assert stats.summary == {
            "alignment_length": 6,
            "output_length": 3,
            "trimmed_length": 3,
            "trimmed_percentage": 50.0,
        }
        assert isinstance(trim_run.version, str)

    def test_threads_must_be_positive(self):
        with pytest.raises(ValueError, match="threads must be an integer >= 1"):
            clipkit(
                input_file_path="tests/integration/samples/simple.fa",
                mode=TrimmingMode.gappy,
                gaps=0.3,
                sequence_type="nt",
                threads=0,
            )

    def test_requires_exactly_one_input_source(self):
        with pytest.raises(
            ValueError,
            match="Provide exactly one of raw_alignment or input_file_path.",
        ):
            clipkit(
                raw_alignment=">1\nA\n>2\nA\n",
                input_file_path="tests/integration/samples/simple.fa",
            )

        with pytest.raises(
            ValueError,
            match="Provide exactly one of raw_alignment or input_file_path.",
        ):
            clipkit()

    def test_empty_raw_alignment_rejected(self):
        with pytest.raises(ValueError, match="raw_alignment cannot be empty."):
            clipkit(raw_alignment="")

    def test_empty_input_file_path_rejected(self):
        with pytest.raises(ValueError, match="input_file_path cannot be empty."):
            clipkit(input_file_path="")

    def test_invalid_sequence_type_rejected(self):
        with pytest.raises(ValueError, match="sequence_type must be one of"):
            clipkit(
                input_file_path="tests/integration/samples/simple.fa",
                sequence_type="protein",
            )

    def test_plot_trim_report_path_writes_html(self, tmp_path):
        report_path = tmp_path / "api_trim_report.html"
        trim_run, stats = clipkit(
            input_file_path="tests/integration/samples/simple.fa",
            mode=TrimmingMode.gappy,
            gaps=0.3,
            sequence_type="nt",
            plot_trim_report_path=str(report_path),
        )
        assert stats.summary["trimmed_length"] == 2
        assert isinstance(trim_run.trimmed, MultipleSeqAlignment)
        assert report_path.exists()

    def test_heterotachy_mode(self):
        trim_run, stats = clipkit(
            input_file_path="tests/integration/samples/simple.fa",
            mode=TrimmingMode.heterotachy,
            gaps=0.8,
            sequence_type="nt",
        )

        assert stats.summary == {
            "alignment_length": 6,
            "output_length": 4,
            "trimmed_length": 2,
            "trimmed_percentage": 33.333,
        }
        assert isinstance(trim_run.version, str)

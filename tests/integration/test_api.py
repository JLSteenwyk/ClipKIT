import pytest

from Bio.Align import MultipleSeqAlignment
from clipkit import clipkit
from clipkit.files import FileFormat
from clipkit.modes import StopCodonMode, TrimmingMode
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

    @pytest.mark.parametrize(
        "remove_stop_codons, expected_first_sequence, terminal, internal",
        [
            (StopCodonMode.terminal, "ATGTAAACC---", 1, 0),
            ("internal", "ATG---ACCTGA", 0, 1),
            (StopCodonMode.all, "ATG---ACC---", 1, 1),
        ],
    )
    def test_stop_codon_masking_api(
        self,
        remove_stop_codons,
        expected_first_sequence,
        terminal,
        internal,
    ):
        trim_run, _ = clipkit(
            raw_alignment=">dna\nATGTAAACCTGA\n>control\nATGCAAACCCAA\n",
            mode=TrimmingMode.gappy,
            gaps=0.9,
            codon=True,
            sequence_type="nt",
            remove_stop_codons=remove_stop_codons,
        )

        assert str(trim_run.trimmed[0].seq) == expected_first_sequence
        assert trim_run.stop_codon_masking.terminal_masked == terminal
        assert trim_run.stop_codon_masking.internal_masked == internal

    def test_stop_codon_masking_precedes_codon_gap_trimming(self):
        trim_run, stats = clipkit(
            raw_alignment=">stop\nATGTAA\n>control\nATGCAA\n",
            mode=TrimmingMode.gappy,
            gaps=0.5,
            codon=True,
            sequence_type="nt",
            remove_stop_codons="terminal",
        )

        assert [str(record.seq) for record in trim_run.trimmed] == ["ATG", "ATG"]
        assert stats.output_length == 3
        assert trim_run.stop_codon_masking.terminal_masked == 1

    def test_omitting_stop_codon_option_preserves_existing_behavior(self):
        trim_run, _ = clipkit(
            raw_alignment=">stop\nATGTAA\n>control\nATGCAA\n",
            mode=TrimmingMode.gappy,
            gaps=0.9,
            codon=True,
            sequence_type="nt",
        )

        assert [str(record.seq) for record in trim_run.trimmed] == [
            "ATGTAA",
            "ATGCAA",
        ]
        assert trim_run.stop_codon_masking.summary == {
            "mode": None,
            "terminal_masked": 0,
            "internal_masked": 0,
            "total_masked": 0,
        }

    @pytest.mark.parametrize(
        "kwargs, message",
        [
            ({"remove_stop_codons": "terminal"}, "codon-aware trimming"),
            (
                {
                    "codon": True,
                    "sequence_type": "aa",
                    "remove_stop_codons": "terminal",
                },
                "nucleotide input",
            ),
            (
                {
                    "raw_alignment": ">one\nATGTA\n>two\nATGCA\n",
                    "codon": True,
                    "remove_stop_codons": "terminal",
                },
                "alignment length divisible by 3",
            ),
            (
                {"codon": True, "remove_stop_codons": "unsupported"},
                "remove_stop_codons must be one of",
            ),
        ],
    )
    def test_stop_codon_masking_validation(self, kwargs, message):
        options = {
            "raw_alignment": ">one\nATGTAA\n>two\nATGCAA\n",
            "mode": TrimmingMode.gappy,
            "gaps": 0.9,
            "sequence_type": "nt",
        }
        options.update(kwargs)

        with pytest.raises(ValueError, match=message):
            clipkit(**options)

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

    def test_output_file_path_writes_alignment_and_returns_path(self, tmp_path):
        output_path = tmp_path / "api_trimmed.fa"

        returned_path, stats = clipkit(
            input_file_path="tests/integration/samples/simple.fa",
            output_file_path=str(output_path),
            output_file_format=FileFormat.fasta,
            mode=TrimmingMode.gappy,
            gaps=0.3,
            sequence_type="nt",
        )

        assert returned_path == str(output_path)
        assert stats.summary["trimmed_length"] == 2
        assert output_path.read_text() == (
            ">1\nAGAT\n>2\nAGAT\n>3\nAGTA\n>4\nAATA\n>5\nAaT-\n"
        )

    def test_complement_property_exposes_trimmed_columns(self):
        trim_run, _ = clipkit(
            input_file_path="tests/integration/samples/simple.fa",
            mode=TrimmingMode.gappy,
            gaps=0.3,
            sequence_type="nt",
        )

        assert [str(record.seq) for record in trim_run.complement] == [
            "-T",
            "--",
            "--",
            "G-",
            "C-",
        ]

    def test_entropy_mode_uses_documented_default_threshold(self):
        trim_run, stats = clipkit(
            input_file_path="tests/integration/samples/simple.fa",
            mode=TrimmingMode.entropy,
            gaps=None,
            sequence_type="nt",
        )

        assert trim_run.gaps == 0.8
        assert stats.summary == {
            "alignment_length": 6,
            "output_length": 2,
            "trimmed_length": 4,
            "trimmed_percentage": 66.667,
        }

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

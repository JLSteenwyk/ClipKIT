import json

from clipkit.api import clipkit
from clipkit.clipkit import execute
from clipkit.files import FileFormat
from clipkit.logger import logger, log_file_logger
from clipkit.modes import TrimmingMode


def test_execute_invalid_input_format_does_not_crash(tmp_path):
    input_file = tmp_path / "not_alignment.txt"
    input_file.write_text("hello world\n")
    output_file = tmp_path / "out.fa"

    execute(
        input_file=str(input_file),
        input_file_format=None,
        output_file=str(output_file),
        output_file_format=FileFormat.fasta,
        sequence_type=None,
        gaps=0.9,
        gap_characters=None,
        complement=False,
        codon=False,
        ends_only=False,
        mode=TrimmingMode.gappy,
        use_log=False,
        quiet=True,
    )

    assert not output_file.exists()


def test_execute_restores_logger_disabled_flag(tmp_path):
    output_file = tmp_path / "simple.out.fa"
    original_disabled = logger.disabled
    logger.disabled = False

    execute(
        input_file="tests/integration/samples/simple.fa",
        input_file_format=None,
        output_file=str(output_file),
        output_file_format=FileFormat.fasta,
        sequence_type=None,
        gaps=0.9,
        gap_characters=None,
        complement=False,
        codon=False,
        ends_only=False,
        mode=TrimmingMode.gappy,
        use_log=False,
        quiet=True,
    )

    assert logger.disabled is False
    logger.disabled = original_disabled


def test_execute_cleans_up_log_file_handler(tmp_path):
    output_file = tmp_path / "simple.logtest.fa"
    initial_handlers = len(log_file_logger.handlers)

    execute(
        input_file="tests/integration/samples/simple.fa",
        input_file_format=None,
        output_file=str(output_file),
        output_file_format=FileFormat.fasta,
        sequence_type=None,
        gaps=0.9,
        gap_characters=None,
        complement=False,
        codon=False,
        ends_only=False,
        mode=TrimmingMode.gappy,
        use_log=True,
        quiet=True,
    )

    assert len(log_file_logger.handlers) == initial_handlers


def test_api_restores_logger_disabled_flag():
    original_disabled = logger.disabled
    logger.disabled = False

    clipkit(
        input_file_path="tests/integration/samples/simple.fa",
        mode=TrimmingMode.gappy,
        gaps=0.3,
        sequence_type="nt",
    )

    assert logger.disabled is False
    logger.disabled = original_disabled


def test_execute_dry_run_writes_no_output_files(tmp_path):
    output_file = tmp_path / "dryrun.out.fa"

    execute(
        input_file="tests/integration/samples/simple.fa",
        input_file_format=None,
        output_file=str(output_file),
        output_file_format=FileFormat.fasta,
        sequence_type=None,
        gaps=0.9,
        gap_characters=None,
        complement=True,
        codon=False,
        ends_only=False,
        mode=TrimmingMode.gappy,
        use_log=True,
        quiet=True,
        dry_run=True,
    )

    assert not output_file.exists()
    assert not (tmp_path / "dryrun.out.fa.complement").exists()
    assert not (tmp_path / "dryrun.out.fa.log").exists()


def test_execute_validate_only_writes_no_output_files(tmp_path):
    output_file = tmp_path / "validate.out.fa"

    execute(
        input_file="tests/integration/samples/simple.fa",
        input_file_format=None,
        output_file=str(output_file),
        output_file_format=FileFormat.fasta,
        sequence_type=None,
        gaps=0.9,
        gap_characters=None,
        complement=True,
        codon=False,
        ends_only=False,
        mode=TrimmingMode.gappy,
        use_log=True,
        quiet=True,
        validate_only=True,
    )

    assert not output_file.exists()
    assert not (tmp_path / "validate.out.fa.complement").exists()
    assert not (tmp_path / "validate.out.fa.log").exists()


def test_execute_writes_report_json_for_dry_run(tmp_path):
    output_file = tmp_path / "report.out.fa"
    report_file = tmp_path / "report.json"

    execute(
        input_file="tests/integration/samples/simple.fa",
        input_file_format=None,
        output_file=str(output_file),
        output_file_format=FileFormat.fasta,
        sequence_type=None,
        gaps=0.9,
        gap_characters=None,
        complement=False,
        codon=False,
        ends_only=False,
        mode=TrimmingMode.gappy,
        use_log=False,
        quiet=True,
        dry_run=True,
        report_json=str(report_file),
    )

    assert report_file.exists()
    payload = json.loads(report_file.read_text())
    assert payload["status"] == "completed"
    assert payload["dry_run"] is True
    assert payload["validate_only"] is False
    assert "stats" in payload


def test_execute_writes_report_json_for_validate_only(tmp_path):
    output_file = tmp_path / "validate.report.out.fa"
    report_file = tmp_path / "validate.report.json"

    execute(
        input_file="tests/integration/samples/simple.fa",
        input_file_format=None,
        output_file=str(output_file),
        output_file_format=FileFormat.fasta,
        sequence_type=None,
        gaps=0.9,
        gap_characters=None,
        complement=False,
        codon=False,
        ends_only=False,
        mode=TrimmingMode.gappy,
        use_log=False,
        quiet=True,
        validate_only=True,
        report_json=str(report_file),
    )

    assert report_file.exists()
    payload = json.loads(report_file.read_text())
    assert payload["status"] == "validated"
    assert payload["validate_only"] is True
    assert payload["dry_run"] is False
    assert "stats" not in payload

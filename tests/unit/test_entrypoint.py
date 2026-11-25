import subprocess
import sys
from pathlib import Path

import pytest

from clipkit.clipkit import main


def _run_clipkit(*args, **kwargs):
    cmd = [sys.executable, "-m", "clipkit"] + list(args)
    return subprocess.run(cmd, **kwargs)


class TestEntrypoint(object):
    def test_help(self):
        result = _run_clipkit("--help")
        assert result.returncode == 0

    def test_run(self):
        sample = Path("tests/integration/samples/simple.fa")
        result = _run_clipkit(str(sample))
        assert result.returncode == 0

    def test_input_error(self):
        result = _run_clipkit("/file/doesnt/exist", capture_output=True, text=True)
        assert result.returncode == 0
        assert "Input file does not exist" in result.stdout

    def test_run_no_args(self):
        result = _run_clipkit()
        assert result.returncode == 0

    def test_remove_stop_codons_cli(self, tmp_path):
        input_file = tmp_path / "stops.fa"
        output_file = tmp_path / "stops.out.fa"
        input_file.write_text(">stop\nATGTAATGA\n>control\nATGCAACAA\n")

        result = _run_clipkit(
            str(input_file),
            "--output",
            str(output_file),
            "--mode",
            "gappy",
            "--gaps",
            "0.9",
            "--sequence_type",
            "nt",
            "--codon",
            "--remove_stop_codons",
            "all",
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert output_file.read_text() == ">stop\nATG------\n>control\nATGCAACAA\n"
        assert "Stop codon masking mode: all" in result.stdout
        assert "Terminal stop codons masked: 1" in result.stdout
        assert "Internal stop codons masked: 1" in result.stdout

    def test_remove_stop_codons_cli_requires_codon_flag(self, tmp_path):
        input_file = tmp_path / "stops.fa"
        input_file.write_text(">stop\nATGTAA\n>control\nATGCAA\n")

        result = _run_clipkit(
            str(input_file),
            "--remove_stop_codons",
            "terminal",
            capture_output=True,
            text=True,
        )

        assert result.returncode == 2
        assert "Stop codon masking requires --codon" in result.stderr

    def test_main_without_arguments_displays_help(self, capsys):
        main([])

        captured = capsys.readouterr()
        assert "Usage: clipkit <input> [optional arguments]" in captured.err

    def test_main_reports_stop_codon_validation_errors(self, tmp_path, capsys):
        input_file = tmp_path / "stops.fa"
        input_file.write_text(">stop\nATGTAA\n>control\nATGCAA\n")

        with pytest.raises(SystemExit) as error:
            main([str(input_file), "--remove_stop_codons", "terminal"])

        assert error.value.code == 2
        assert "Stop codon masking requires --codon" in capsys.readouterr().err


def test_information_flags_leave_numeric_runtime_unloaded():
    script = """
import sys
from clipkit.cli import main

for flag in ("--help", "--version"):
    try:
        main([flag])
    except SystemExit as error:
        assert error.code == 0
assert "numpy" not in sys.modules
assert "Bio" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", script], check=True, capture_output=True)


def test_runtime_import_preserves_api_callable_and_enum_pickle_paths():
    script = """
import importlib
import pickle
import clipkit as package

assert "clipkit" in dir(package)
engine = importlib.import_module("clipkit.clipkit")
assert callable(engine.execute)

from clipkit import clipkit
from clipkit.api import clipkit as public_api
from clipkit.files import FileFormat
from clipkit.helpers import SeqType

assert clipkit is public_api and callable(clipkit)
assert FileFormat.__module__ == "clipkit.files"
assert SeqType.__module__ == "clipkit.helpers"
for option in (FileFormat.fasta, SeqType.nt):
    assert pickle.loads(pickle.dumps(option)) is option
"""
    subprocess.run([sys.executable, "-c", script], check=True, capture_output=True)

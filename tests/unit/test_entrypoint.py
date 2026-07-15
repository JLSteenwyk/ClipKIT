import subprocess
import sys
from pathlib import Path


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

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

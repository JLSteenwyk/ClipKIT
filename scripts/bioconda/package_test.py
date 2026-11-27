"""Smoke-test the installed ClipKIT Conda package."""

from importlib.metadata import version
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from Bio import AlignIO


def main():
    release = version("clipkit")
    assert subprocess.check_output(["clipkit", "--version"], text=True).strip() == f"clipkit {release}"
    with TemporaryDirectory() as directory:
        source = Path(directory) / "input.fasta"
        output = Path(directory) / "output.fasta"
        source.write_text(">first\nAC-T\n>second\nAC-T\n")
        subprocess.run(
            ["clipkit", str(source), "-m", "gappy", "-o", str(output)], check=True
        )
        alignment = AlignIO.read(output, "fasta")
        assert [(record.id, str(record.seq)) for record in alignment] == [
            ("first", "ACT"),
            ("second", "ACT"),
        ]


if __name__ == "__main__":
    main()

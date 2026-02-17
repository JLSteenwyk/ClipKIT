import os

from Bio import AlignIO

from clipkit.msa import MSA
from clipkit.modes import TrimmingMode
from clipkit.plot_report import write_trim_plot_report


def get_biopython_msa(file_path, file_format="fasta"):
    return AlignIO.read(open(file_path), file_format)


def test_write_trim_plot_report(tmp_path):
    bio_msa = get_biopython_msa("tests/unit/examples/simple.fa")
    msa = MSA.from_bio_msa(bio_msa)
    msa.trim(mode=TrimmingMode.gappy, gap_threshold=0.3)

    output_file = tmp_path / "trim_report.html"
    write_trim_plot_report(str(output_file), msa, mode="gappy", gaps=0.3)

    assert os.path.isfile(output_file)
    content = output_file.read_text(encoding="utf-8")
    assert "ClipKIT Trim Report" in content
    assert "trimmed" in content
    assert "payload" in content
    assert "10.1371/journal.pbio.3001007" in content
    assert "journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007" in content
    assert "Export Per-site Tracks (PNG)" in content
    assert "Export Alignment Preview (PNG)" in content

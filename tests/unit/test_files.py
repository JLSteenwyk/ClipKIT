import pytest
from pathlib import Path

from Bio import AlignIO
from clipkit.files import get_alignment_and_format, FileFormat
from clipkit.files import get_custom_sites_to_trim

here = Path(__file__)
SAMPLES_DIR = here.parent.parent / "integration" / "samples"


class TestAutomaticFileTypeDetermination(object):
    def test_get_alignment_and_format_when_format_is_provided(self):
        # setup
        in_file = f"{here.parent}/examples/simple.fa"
        file_format = "fasta"

        # execution
        alignment, in_file_format = get_alignment_and_format(in_file, file_format)

        # check results
        assert in_file_format == FileFormat.fasta
        assert alignment.get_alignment_length() == 6

    def test_get_alignment_and_format_when_format_is_not_provided(self):
        # setup
        in_file = f"{here.parent}/examples/simple.fa"
        file_format = None

        # execution
        alignment, in_file_format = get_alignment_and_format(in_file, file_format)

        # check results
        assert in_file_format == FileFormat.fasta
        assert alignment.get_alignment_length() == 6

    def test_get_alignment_and_format_raises_error_when_file_not_known(self, mocker):
        in_file = ""
        file_format = None
        mocker.patch("clipkit.files.AlignIO.read", side_effect=ValueError())

        with pytest.raises(Exception) as excinfo:
            get_alignment_and_format(in_file, file_format)
        assert "No such file or directory" in str(excinfo.value)

    def test_get_alignment_and_format_raises_error_when_detection_fails(self, mocker):
        in_file = f"{here.parent}/examples/simple.fa"
        file_format = None
        mocker.patch("clipkit.files.AlignIO.read", side_effect=ValueError())

        with pytest.raises(Exception) as excinfo:
            get_alignment_and_format(in_file, file_format)
        assert "File could not be read" in str(excinfo.value)

    def test_get_alignment_and_format_detects_ecomp(self):
        in_file = SAMPLES_DIR / "12_YIL115C_Anc_2.253_aa_aln.ecomp"

        alignment, in_file_format = get_alignment_and_format(str(in_file), None)

        assert in_file_format == FileFormat.ecomp
        assert alignment.get_alignment_length() > 0
        assert alignment.annotations.get("ecomp_metadata")

    def test_get_alignment_and_format_supports_explicit_ecomp(self):
        in_file = SAMPLES_DIR / "12_YIL115C_Anc_2.253_aa_aln.ecomp"

        alignment, in_file_format = get_alignment_and_format(
            str(in_file), FileFormat.ecomp.value
        )

        assert in_file_format == FileFormat.ecomp
        assert alignment.get_alignment_length() > 0


class TestCustomSitesParsing(object):
    def test_get_custom_sites_to_trim_invalid_format(self, tmp_path):
        cst_file = tmp_path / "bad.cst"
        cst_file.write_text("1\n")

        with pytest.raises(ValueError, match="Invalid CST format"):
            get_custom_sites_to_trim(str(cst_file), aln_length=10)

    def test_get_custom_sites_to_trim_out_of_range(self, tmp_path):
        cst_file = tmp_path / "bad.cst"
        cst_file.write_text("11\ttrim\n")

        with pytest.raises(ValueError, match="out of range"):
            get_custom_sites_to_trim(str(cst_file), aln_length=10)

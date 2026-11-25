import pytest
import numpy as np
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import clipkit.files as files_module
from clipkit.ecomp.reader import HEADER_MAGIC
from clipkit.exceptions import InvalidInputFileFormat
from clipkit.files import get_alignment_and_format, FileFormat
from clipkit.files import get_custom_sites_to_trim
from clipkit.helpers import write_msa, write_complement
from clipkit.msa import MSA

here = Path(__file__)
SAMPLES_DIR = here.parent.parent / "integration" / "samples"


@pytest.mark.parametrize(
    ("data", "falls_back"),
    [
        (b">one first\nAC GT\nTA\n>two second\nGT\nACAA", False),
        (b">one first\r\nACGTTA\r\n>two second\r\nGTACAA\r\n", True),
        (b">one > first\nACGTTA\n>two second\nGTACAA", False),
        (b">one first\nAC>GTTA\n>two second\nGTACAAA", False),
    ],
)
def test_block_fasta_preserves_biopython_records(tmp_path, monkeypatch, mocker, data, falls_back):
    path = tmp_path / "alignment.fa"
    path.write_bytes(data)
    with path.open() as handle:
        expected = AlignIO.read(handle, "fasta")
    monkeypatch.setattr(files_module, "_FASTA_BLOCK_BYTES", 20)
    monkeypatch.setattr(files_module, "_FASTA_BLOCK_MIN_FILE_BYTES", 20)
    biopython_read = mocker.spy(AlignIO, "read")

    actual, format_detected = get_alignment_and_format(path, None)

    assert format_detected == FileFormat.fasta
    assert biopython_read.called == falls_back
    assert [(r.seq, r.id, r.name, r.description) for r in actual] == [
        (r.seq, r.id, r.name, r.description) for r in expected
    ]


@pytest.mark.parametrize("repeats", [300, 2101])
def test_batched_fasta_matches_biopython_for_kept_and_complementary_sites(tmp_path, repeats):
    alignment = MultipleSeqAlignment([
        SeqRecord(Seq("AaX-" * repeats), id="one", description="one full header"),
        SeqRecord(Seq("gTc?" * repeats), id="two", description="two\nheader"),
    ])
    msa = MSA.from_bio_msa(alignment, sequence_type="aa")
    msa.trim(site_positions_to_trim=list(range(2 * repeats - 1)))
    path = tmp_path / "actual.fa"
    write_msa(msa, path, FileFormat.fasta)
    write_complement(msa, path, FileFormat.fasta)

    for result, records in (
        (path, msa.to_bio_msa()),
        (path.with_suffix(".fa.complement"), msa.complement_to_bio_msa()),
    ):
        expected = tmp_path / "expected.fa"
        SeqIO.write(records, expected, "fasta")
        assert result.read_bytes() == expected.read_bytes()


def test_wide_fasta_validation_preserves_destination_on_bad_sequence(tmp_path):
    records = np.full((2, 4097), "A", dtype="U1")
    records[1, -1] = "é"
    msa = MSA([{"description": "first"}, {"description": "second"}], records)
    output = tmp_path / "existing.fa"
    output.write_text("preserve this file\n")

    with pytest.raises(UnicodeEncodeError):
        write_msa(msa, output, FileFormat.fasta)

    assert output.read_text() == "preserve this file\n"


def test_wide_fasta_validates_lengths_after_removing_null_padding(tmp_path):
    records = np.full((2, 4097), "A", dtype="U1")
    records[1, -1] = ""
    msa = MSA([{"description": "first"}, {"description": "second"}], records)
    output = tmp_path / "existing.fa"
    output.write_text("preserve this file\n")

    with pytest.raises(ValueError, match="same length"):
        write_msa(msa, output, FileFormat.fasta)

    assert output.read_text() == "preserve this file\n"


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

    def test_corrupt_ecomp_header_falls_through_to_format_detection(self, tmp_path):
        corrupt_archive = tmp_path / "corrupt.ecomp"
        corrupt_archive.write_bytes(HEADER_MAGIC + b"not-an-ecomp-archive")

        with pytest.raises(InvalidInputFileFormat, match="File could not be read"):
            get_alignment_and_format(str(corrupt_archive), None)

    def test_format_detection_skips_assertion_errors(self, mocker):
        input_file = f"{here.parent}/examples/simple.fa"
        mocker.patch("clipkit.files.AlignIO.read", side_effect=AssertionError())

        with pytest.raises(InvalidInputFileFormat, match="File could not be read"):
            get_alignment_and_format(input_file, None)


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

    def test_get_custom_sites_to_trim_keep_only_ignores_blank_lines(self, tmp_path):
        cst_file = tmp_path / "keep-only.cst"
        cst_file.write_text("\n4\tkeep\n\n2\tkeep\n4\tkeep\n")

        assert get_custom_sites_to_trim(str(cst_file), aln_length=5) == [0, 2, 4]

    def test_get_custom_sites_to_trim_rejects_non_numeric_position(self, tmp_path):
        cst_file = tmp_path / "bad-position.cst"
        cst_file.write_text("first\ttrim\n")

        with pytest.raises(ValueError, match="Invalid CST position"):
            get_custom_sites_to_trim(str(cst_file), aln_length=5)

    def test_get_custom_sites_to_trim_rejects_unknown_action(self, tmp_path):
        cst_file = tmp_path / "bad-action.cst"
        cst_file.write_text("1\tdrop\n")

        with pytest.raises(ValueError, match="Invalid CST action"):
            get_custom_sites_to_trim(str(cst_file), aln_length=5)

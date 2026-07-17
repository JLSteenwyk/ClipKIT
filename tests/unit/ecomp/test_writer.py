import struct

from Bio import AlignIO

from clipkit.ecomp import read_ecomp, write_ecomp
from clipkit.ecomp.reader import (
    HEADER_MAGIC,
    HEADER_STRUCT,
    INLINE_METADATA_VERSION,
)


def test_write_ecomp_sanitizes_inherited_metadata_and_roundtrips(tmp_path):
    alignment = AlignIO.read("tests/unit/examples/simple.fa", "fasta")
    output_path = tmp_path / "alignment.ecomp"
    inherited_metadata = {
        "format_version": (0, 2, 0),
        "source_format": "fasta",
        "sequence_permutation": [1, 0],
        "run_length_blocks": [[0, 2]],
        "payload_encoding": "raw",
        "fallback": {"type": "raw"},
    }

    metadata = write_ecomp(alignment, output_path, inherited_metadata)
    decoded_alignment, decoded_metadata = read_ecomp(output_path)

    assert metadata["format_version"] == "0.2.0"
    assert "sequence_permutation" not in metadata
    assert "run_length_blocks" not in metadata
    assert metadata["payload_encoding"] == "gzip"
    assert decoded_metadata == metadata
    assert [str(record.seq) for record in decoded_alignment] == [
        str(record.seq) for record in alignment
    ]


def test_write_ecomp_falls_back_to_supported_header_for_bad_version(tmp_path):
    alignment = AlignIO.read("tests/unit/examples/simple.fa", "fasta")
    output_path = tmp_path / "bad-version.ecomp"

    write_ecomp(alignment, output_path, {"format_version": "not.a.version"})

    header_size = struct.calcsize(HEADER_STRUCT)
    header = struct.unpack(HEADER_STRUCT, output_path.read_bytes()[:header_size])
    assert header[0] == HEADER_MAGIC
    assert header[1:4] == INLINE_METADATA_VERSION

import base64
import gzip
import json
import lzma
import struct
import types
import zlib
from hashlib import sha256
from pathlib import Path

import pytest

from clipkit.ecomp.reader import (
    HEADER_MAGIC,
    LEGACY_HEADER_STRUCT,
    METADATA_COMPRESSED_MAGIC,
    METADATA_CODEC_VERSION,
    METADATA_SUFFIX,
    PERM_MAGIC,
    PERM_VERSION,
    SEQ_ID_MAGIC,
    SEQ_ID_VERSION,
    INLINE_METADATA_VERSION,
    EcompDecodeError,
    RunLengthBlock,
    decode_blocks,
    read_ecomp,
    _canonical_code_maps,
    _decode_bitmask,
    _decode_metadata,
    _decode_permutation,
    _decode_residue_stream,
    _decode_residues,
    _decode_sequence_ids,
    _alignment_from_fasta_bytes,
    _decode_varint,
    _decompress_alignment,
    _derive_metadata_path,
    _extract_permutation_chunk,
    _iter_deviation_indices,
    _pack_codes,
    _popcount,
    _read_archive,
    _read_metadata,
    _read_varint,
    _unpack_codes,
)


def encode_varint(value: int) -> bytes:
    parts = []
    while True:
        byte = value & 0x7F
        value >>= 7
        if value:
            parts.append(byte | 0x80)
        else:
            parts.append(byte)
            break
    return bytes(parts)


def write_inline_archive(
    path: Path, payload: bytes, metadata: dict[str, object], version=(0, 2, 0)
) -> None:
    header = struct.pack(
        LEGACY_HEADER_STRUCT,
        HEADER_MAGIC,
        version[0],
        version[1],
        version[2],
        len(payload),
    )
    metadata_bytes = json.dumps(metadata).encode("utf-8")
    archive = header + struct.pack(">Q", len(metadata_bytes))
    archive += payload + metadata_bytes
    path.write_bytes(archive)


def write_legacy_archive(
    path: Path, payload: bytes, metadata: dict[str, object], version=(0, 1, 0)
) -> Path:
    header = struct.pack(
        LEGACY_HEADER_STRUCT,
        HEADER_MAGIC,
        version[0],
        version[1],
        version[2],
        len(payload),
    )
    path.write_bytes(header + payload)
    meta_path = path.with_suffix(METADATA_SUFFIX)
    meta_path.write_text(json.dumps(metadata))
    return meta_path


def build_seq_id_block(ids: list[str], mode: int = 0) -> bytes:
    inner = encode_varint(len(ids))
    for name in ids:
        encoded = name.encode("utf-8")
        inner += encode_varint(len(encoded)) + encoded
    if mode == 0:
        payload = inner
    elif mode == 2:
        payload = zlib.compress(inner)
    else:
        payload = inner
    block = bytes([mode]) + payload
    header = SEQ_ID_MAGIC + bytes([SEQ_ID_VERSION])
    header += encode_varint(len(block))
    return header + block


def build_permutation_payload(
    permutation: list[int], compressed: bool = False
) -> tuple[bytes, dict[str, object]]:
    width = 1
    payload = b"".join(int(value).to_bytes(width, "little") for value in permutation)
    payload_bytes = zlib.compress(payload) if compressed else payload
    chunk = bytearray(PERM_MAGIC)
    chunk.append(PERM_VERSION)
    flags = 0x01 if compressed else 0
    chunk.append(flags)
    chunk.extend(encode_varint(len(permutation)))
    chunk.extend(encode_varint(len(payload_bytes)))
    chunk.extend(payload_bytes)
    meta = {"encoding": "payload", "length": len(chunk)}
    return bytes(chunk), meta


class TestDecodeMetadata:
    def test_compressed_metadata_roundtrip(self):
        payload = json.dumps({"foo": "bar"}).encode("utf-8")
        compressed = (
            METADATA_COMPRESSED_MAGIC
            + bytes([METADATA_CODEC_VERSION])
            + zlib.compress(payload)
        )
        decoded = _decode_metadata(compressed)
        assert decoded == {"foo": "bar"}

    def test_compressed_metadata_header_truncated(self):
        with pytest.raises(EcompDecodeError):
            _decode_metadata(METADATA_COMPRESSED_MAGIC)

    def test_unsupported_compressed_version(self):
        payload = (
            METADATA_COMPRESSED_MAGIC
            + bytes([METADATA_CODEC_VERSION + 1])
            + b"payload"
        )
        with pytest.raises(EcompDecodeError):
            _decode_metadata(payload)

    def test_malformed_json_raises(self):
        with pytest.raises(EcompDecodeError):
            _decode_metadata(b"{not json}")

    def test_plain_metadata_json(self):
        payload = json.dumps({"baz": 3}).encode("utf-8")
        assert _decode_metadata(payload) == {"baz": 3}


class TestMetadataHelpers:
    def test_read_metadata_from_file(self, tmp_path):
        path = tmp_path / "meta.json"
        payload = json.dumps({"spam": 1}).encode("utf-8")
        path.write_bytes(payload)
        assert _read_metadata(path) == {"spam": 1}

    def test_derive_metadata_path(self):
        base = Path("archive.ecomp")
        derived = _derive_metadata_path(base)
        assert derived.name == "archive.json"


class TestDecodeBitmask:
    def test_no_deviation_returns_zero_mask(self):
        result = _decode_bitmask(0, b"\xFF", 0, 2)
        assert result == b"\x00\x00"

    def test_copy_mode_preserves_payload(self):
        result = _decode_bitmask(0, b"\x0F", 3, 2)
        assert result == b"\x0F\x00"

    def test_position_mode_expands_positions(self):
        positions = [0, 3, 9]
        encoded_positions = b"".join(
            encode_varint(pos - prev)
            for prev, pos in zip([-1] + positions[:-1], positions)
        )
        payload = encode_varint(len(encoded_positions)) + encoded_positions
        result = _decode_bitmask(1, payload, len(positions), 2)
        assert result == bytes([0b00001001, 0b00000010])

    def test_position_mode_rejects_out_of_range(self):
        encoded_positions = encode_varint(17)
        payload = encode_varint(len(encoded_positions)) + encoded_positions
        with pytest.raises(EcompDecodeError):
            _decode_bitmask(1, payload, 1, 2)

    def test_rle_mode_rehydrates_sequence(self):
        payload = bytes([0b00000111, 2]) + bytes([0, 1])
        result = _decode_bitmask(2, payload, 4, 3)
        assert result == bytes([0b00000111, 0b00000111, 0])

    def test_decode_bitmask_unknown_mode(self):
        with pytest.raises(EcompDecodeError):
            _decode_bitmask(9, b"", 1, 1)

    def test_decode_bitmask_rle_truncated(self):
        with pytest.raises(EcompDecodeError):
            _decode_bitmask(2, bytes([0xAA]), 1, 1)

    def test_decode_bitmask_rle_pads_with_zeros(self):
        payload = bytes([0b00000001, 1])
        result = _decode_bitmask(2, payload, 1, 3)
        assert result.endswith(b"\x00\x00")


class TestArchiveReading:
    def test_read_ecomp_uses_fallback(self, tmp_path):
        fasta = b"\n".join([b">s0", b"AC", b">s1", b"AG"]) + b"\n"
        payload = gzip.compress(fasta)
        metadata = {
            "fallback": {"type": "gzip"},
            "source_format": "fasta",
        }
        archive = tmp_path / "fallback.ecomp"
        write_inline_archive(archive, payload, metadata)

        alignment, decoded = read_ecomp(archive)

        assert decoded == metadata
        assert alignment.annotations["ecomp_metadata"] == metadata
        assert alignment.annotations["ecomp_version"] == INLINE_METADATA_VERSION
        assert [record.id for record in alignment] == ["s0", "s1"]
        assert alignment.get_alignment_length() == 2

    def test_read_archive_inline_metadata(self, tmp_path):
        payload = b"data"
        metadata = {"meta": 1}
        archive = tmp_path / "inline.ecomp"
        write_inline_archive(archive, payload, metadata)

        extracted_payload, decoded_meta, version = _read_archive(archive)

        assert extracted_payload == payload
        assert decoded_meta == metadata
        assert version == INLINE_METADATA_VERSION

    def test_read_archive_legacy_requires_sidecar(self, tmp_path):
        payload = b""
        archive = tmp_path / "legacy.ecomp"
        header = struct.pack(
            LEGACY_HEADER_STRUCT, HEADER_MAGIC, 0, 1, 0, len(payload)
        )
        archive.write_bytes(header + payload)

        with pytest.raises(FileNotFoundError):
            _read_archive(archive)


class TestDecompressAlignment:
    def test_decompress_alignment_reconstructs_sequences(self, monkeypatch):
        alphabet = ["A", "C", "G"]
        metadata = {
            "alignment_length": 3,
            "num_sequences": 2,
            "sequence_ids": ["s0", "s1"],
            "alphabet": alphabet,
            "bitmask_bytes": 1,
            "payload_encoding": "zlib",
            "sequence_permutation": [1, 0],
        }
        payload = zlib.compress(b"blocks")
        residues = _pack_codes([1], 2)
        blocks = [
            RunLengthBlock("A", bytes([0]), b"", 1),
            RunLengthBlock("A", bytes([0b00000001]), residues, 1),
            RunLengthBlock("G", bytes([0]), b"", 1),
        ]

        def fake_decode(data, bitmask_bytes, bits_per_symbol, alpha):
            assert data == b"blocks"
            assert bitmask_bytes == 1
            assert bits_per_symbol == 2
            assert list(alpha) == alphabet
            return blocks

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        digest = sha256()
        for seq in ("AAG", "ACG"):
            digest.update(seq.encode("utf-8"))
            digest.update(b"\n")
        metadata["checksum_sha256"] = digest.hexdigest()

        alignment = _decompress_alignment(payload, metadata)

        assert [record.id for record in alignment] == ["s1", "s0"]
        assert [str(record.seq) for record in alignment] == ["AAG", "ACG"]
        assert metadata["sequence_permutation"] == [1, 0]

    def test_decompress_alignment_invalid_bitmask(self):
        metadata = {
            "alignment_length": 1,
            "num_sequences": 1,
            "sequence_ids": ["s0"],
            "alphabet": ["A"],
            "bitmask_bytes": 0,
            "bits_per_symbol": 1,
        }
        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_sequence_id_count_mismatch(self):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 2,
            "sequence_ids": ["s0"],
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
        }
        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_missing_bits_per_symbol(self):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 0,
            "sequence_ids": [],
            "alphabet": [],
            "bitmask_bytes": 1,
        }
        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_infers_sequence_ids(self, monkeypatch):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 2,
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
        }
        payload = build_seq_id_block(["s0", "s1"]) + b""

        def fake_decode(data, *_):
            assert data == b""
            return []

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        alignment = _decompress_alignment(payload, metadata)

        assert [record.id for record in alignment] == ["s0", "s1"]
        assert metadata["sequence_ids"] == ["s0", "s1"]

    def test_decompress_alignment_mismatched_sequence_ids(self, monkeypatch):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 2,
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
            "sequence_ids": ["x", "y"],
        }
        payload = build_seq_id_block(["s0", "s1"])

        def fake_decode(data, *_):
            return []

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        with pytest.raises(EcompDecodeError):
            _decompress_alignment(payload, metadata)

    def test_decompress_alignment_missing_sequence_ids(self, monkeypatch):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 1,
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
        }

        def fake_decode(data, *_):
            assert data == b""
            return []

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_columns_exceed_length(self, monkeypatch):
        metadata = {
            "alignment_length": 1,
            "num_sequences": 1,
            "sequence_ids": ["s0"],
            "alphabet": ["A"],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
        }
        blocks = [RunLengthBlock("A", bytes([0]), b"", 2)]

        def fake_decode(*_):
            return blocks

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_columns_too_few(self, monkeypatch):
        metadata = {
            "alignment_length": 2,
            "num_sequences": 1,
            "sequence_ids": ["s0"],
            "alphabet": ["A"],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
        }
        blocks = [RunLengthBlock("A", bytes([0]), b"", 1)]

        def fake_decode(*_):
            return blocks

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_checksum_mismatch(self, monkeypatch):
        metadata = {
            "alignment_length": 1,
            "num_sequences": 1,
            "sequence_ids": ["s0"],
            "alphabet": ["A", "G"],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
            "checksum_sha256": "deadbeef",
        }
        blocks = [RunLengthBlock("A", bytes([0]), b"", 1)]

        def fake_decode(*_):
            return blocks

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_unsupported_encoding(self):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 0,
            "sequence_ids": [],
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
            "payload_encoding": "unknown",
        }
        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_xz_encoding(self, monkeypatch):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 0,
            "sequence_ids": [],
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
            "payload_encoding": "xz",
        }
        payload = lzma.compress(b"rest")

        def fake_decode(data, *_):
            assert data == b"rest"
            return []

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        alignment = _decompress_alignment(payload, metadata)
        assert len(alignment) == 0

    def test_decompress_alignment_raw_encoding(self, monkeypatch):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 0,
            "sequence_ids": [],
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
        }

        def fake_decode(data, *_):
            assert data == b"raw"
            return []

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        alignment = _decompress_alignment(b"raw", metadata)
        assert len(alignment) == 0

    def test_decompress_alignment_requires_zstd(self, monkeypatch):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 0,
            "sequence_ids": [],
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
            "payload_encoding": "zstd",
        }
        monkeypatch.setattr("clipkit.ecomp.reader._ZSTD_AVAILABLE", False)
        monkeypatch.setattr("clipkit.ecomp.reader.zstd", None)
        with pytest.raises(EcompDecodeError):
            _decompress_alignment(b"", metadata)

    def test_decompress_alignment_zstd_path(self, monkeypatch):
        metadata = {
            "alignment_length": 0,
            "num_sequences": 0,
            "sequence_ids": [],
            "alphabet": [],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
            "payload_encoding": "zstd",
        }
        payload = b"data"

        class Dummy:
            def decompress(self, value):
                assert value == payload
                return b"converted"

        dummy_module = types.SimpleNamespace(ZstdDecompressor=lambda: Dummy())
        monkeypatch.setattr("clipkit.ecomp.reader._ZSTD_AVAILABLE", True)
        monkeypatch.setattr("clipkit.ecomp.reader.zstd", dummy_module)

        def fake_decode(data, *_):
            assert data == b"converted"
            return []

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        alignment = _decompress_alignment(payload, metadata)
        assert len(alignment) == 0

    def test_decompress_alignment_permutation_chunk(self, monkeypatch):
        chunk, perm_meta = build_permutation_payload([1, 0], compressed=True)
        metadata = {
            "alignment_length": 1,
            "num_sequences": 2,
            "sequence_ids": ["a", "b"],
            "alphabet": ["A"],
            "bitmask_bytes": 1,
            "bits_per_symbol": 1,
            "sequence_permutation": perm_meta,
        }
        payload = chunk + b"tail"

        def fake_decode(data, *_):
            assert data == b"tail"
            return [RunLengthBlock("A", bytes([0]), b"", 1)]

        monkeypatch.setattr("clipkit.ecomp.reader.decode_blocks", fake_decode)

        alignment = _decompress_alignment(payload, metadata)

        assert metadata["sequence_permutation"] == [1, 0]
        assert [record.id for record in alignment] == ["b", "a"]


    def test_read_archive_legacy_with_sidecar(self, tmp_path):
        metadata = {"foo": "bar"}
        payload = b"abc"
        archive = tmp_path / "legacy_ok.ecomp"
        write_legacy_archive(archive, payload, metadata)

        extracted_payload, decoded_meta, version = _read_archive(archive)

        assert extracted_payload == payload
        assert decoded_meta == metadata
        assert version == (0, 1, 0)

    def test_read_archive_file_too_small(self, tmp_path):
        archive = tmp_path / "small.ecomp"
        archive.write_bytes(b"1234")
        with pytest.raises(EcompDecodeError):
            _read_archive(archive)

    def test_read_archive_invalid_magic(self, tmp_path):
        payload = b""
        archive = tmp_path / "bad_magic.ecomp"
        header = struct.pack(
            LEGACY_HEADER_STRUCT, b"BADMAGIC", 0, 2, 0, len(payload)
        )
        archive.write_bytes(header + payload)
        with pytest.raises(EcompDecodeError):
            _read_archive(archive)

    def test_read_archive_truncated_metadata_header(self, tmp_path):
        payload = b"abc"
        archive = tmp_path / "truncated_meta.ecomp"
        header = struct.pack(
            LEGACY_HEADER_STRUCT, HEADER_MAGIC, 0, 2, 0, len(payload)
        )
        # Declare metadata length but omit bytes
        archive.write_bytes(header + struct.pack(">Q", 5) + payload)
        with pytest.raises(EcompDecodeError):
            _read_archive(archive)

    def test_read_archive_payload_length_exceeds(self, tmp_path):
        payload = b"abc"
        archive = tmp_path / "payload_short.ecomp"
        header = struct.pack(
            LEGACY_HEADER_STRUCT, HEADER_MAGIC, 0, 2, 0, len(payload) + 2
        )
        archive.write_bytes(header + struct.pack(">Q", 0) + payload)
        with pytest.raises(EcompDecodeError):
            _read_archive(archive)

    def test_read_archive_truncated_before_metadata_header(self, tmp_path):
        payload = b""
        archive = tmp_path / "trunc_header.ecomp"
        header = struct.pack(
            LEGACY_HEADER_STRUCT, HEADER_MAGIC, 0, 2, 0, len(payload)
        )
        archive.write_bytes(header)
        with pytest.raises(EcompDecodeError):
            _read_archive(archive)

    def test_alignment_from_fasta_skips_blank_lines(self):
        data = b"\n".join([b">s0", b"AC", b"", b">s1", b"A-"]) + b"\n"
        alignment = _alignment_from_fasta_bytes(data, "fasta")
        assert [record.id for record in alignment] == ["s0", "s1"]

class TestResidueHelpers:
    def test_decode_residue_stream_table_mode(self):
        bitmask = bytes([0b00000101])
        models = {
            "A": {"mode": 0, "residues": ["A", "C", "T"], "width": 2}
        }
        encoded_local = _pack_codes([1, 2], 2)
        alphabet_lookup = {"A": 0, "C": 1, "T": 2, "G": 3}
        result = _decode_residue_stream(
            "A", bitmask, encoded_local, models, 2, alphabet_lookup
        )
        assert result == _pack_codes([1, 2], 2)

    def test_decode_residue_stream_huffman_mode(self):
        residues = ["C", "G"]
        lengths = [1, 2]
        encode_map, decode_map, max_len = _canonical_code_maps(
            residues, lengths
        )
        models = {
            "A": {
                "mode": 1,
                "decode_map": decode_map,
                "max_code_len": max_len,
            }
        }
        bitmask = bytes([0b00000011])
        bits = []
        for residue in ["G", "C"]:
            code, length = encode_map[residue]
            bits.extend(f"{code:0{length}b}")
        while len(bits) % 8:
            bits.append("0")
        encoded = int("".join(bits), 2).to_bytes(len(bits) // 8, "big")
        alphabet_lookup = {"C": 0, "G": 1}
        result = _decode_residue_stream(
            "A", bitmask, encoded, models, 1, alphabet_lookup
        )
        assert result == _pack_codes([1, 0], 1)

    def test_decode_residue_stream_missing_model(self):
        bitmask = bytes([0b00000001])
        with pytest.raises(EcompDecodeError):
            _decode_residue_stream(
                "A", bitmask, b"", {}, 1, {"A": 0}
            )

    def test_decode_residue_stream_missing_alphabet_entry(self):
        bitmask = bytes([0b00000001])
        models = {"A": {"mode": 0, "residues": ["Z"], "width": 1}}
        encoded_local = _pack_codes([0], 1)
        with pytest.raises(EcompDecodeError):
            _decode_residue_stream(
                "A",
                bitmask,
                encoded_local,
                models,
                1,
                {},
            )

    def test_decode_residue_stream_no_deviations(self):
        models = {"A": {"mode": 0, "residues": [], "width": 1}}
        result = _decode_residue_stream(
            "A", bytes([0]), b"", models, 1, {}
        )
        assert result == b""

    def test_decode_residue_stream_table_index_error(self):
        models = {"A": {"mode": 0, "residues": ["C"], "width": 1}}
        encoded_local = _pack_codes([1], 1)
        with pytest.raises(EcompDecodeError):
            _decode_residue_stream(
                "A",
                bytes([0b00000001]),
                encoded_local,
                models,
                1,
                {"C": 0},
            )

    def test_decode_residue_stream_insufficient_huffman_bits(self):
        models = {
            "A": {
                "mode": 1,
                "decode_map": {(1, 0): "C"},
                "max_code_len": 1,
            }
        }
        with pytest.raises(EcompDecodeError):
            _decode_residue_stream(
                "A",
                bytes([0b00000001]),
                b"",
                models,
                1,
                {"C": 0},
            )

    def test_decode_residue_stream_invalid_huffman_code(self):
        models = {
            "A": {
                "mode": 1,
                "decode_map": {},
                "max_code_len": 1,
            }
        }
        with pytest.raises(EcompDecodeError):
            _decode_residue_stream(
                "A",
                bytes([0b00000001]),
                bytes([0x80]),
                models,
                1,
                {"C": 0},
            )

    def test_decode_residues_roundtrip(self):
        data = _pack_codes([0, 2, 1], 2)
        alphabet = ["A", "C", "G"]
        result = _decode_residues(data, 3, 2, alphabet)
        assert result == ["A", "G", "C"]

    def test_decode_residues_insufficient_data(self):
        data = bytes([0b10000000])
        with pytest.raises(EcompDecodeError):
            _decode_residues(data, 2, 6, ["A"])

    def test_decode_residues_exceeds_alphabet(self):
        data = _pack_codes([3], 2)
        with pytest.raises(EcompDecodeError):
            _decode_residues(data, 1, 2, ["A"])


class TestSequenceIdDecoding:
    def make_sequence_buffer(self, mode: int, payload: bytes) -> bytes:
        block = bytes([mode]) + payload
        header = (
            SEQ_ID_MAGIC
            + bytes([SEQ_ID_VERSION])
            + encode_varint(len(block))
        )
        return header + block

    def make_v1_buffer(self, payload: bytes) -> bytes:
        header = SEQ_ID_MAGIC + bytes([1]) + encode_varint(len(payload))
        return header + payload

    def test_decode_sequence_ids_mode_zero(self):
        inner = encode_varint(2)
        for name in (b"seq1", b"seq2"):
            inner += encode_varint(len(name)) + name
        buffer = self.make_sequence_buffer(0, inner)
        ids, remaining = _decode_sequence_ids(buffer + b"tail")
        assert ids == ["seq1", "seq2"]
        assert remaining == b"tail"

    def test_decode_sequence_ids_zlib(self):
        inner = encode_varint(1) + encode_varint(4) + b"name"
        compressed = zlib.compress(inner)
        buffer = self.make_sequence_buffer(2, compressed)
        ids, remaining = _decode_sequence_ids(buffer)
        assert ids == ["name"]
        assert remaining == b""

    def test_decode_sequence_ids_rejects_trailing(self):
        payload = (
            encode_varint(1)
            + encode_varint(4)
            + b"name"
            + b"X"
        )
        buffer = self.make_v1_buffer(payload)
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(buffer)

    def test_decode_sequence_ids_block_truncated(self):
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(SEQ_ID_MAGIC)

    def test_decode_sequence_ids_magic_missing(self):
        header = b"BAD!" + bytes([SEQ_ID_VERSION]) + encode_varint(0)
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(header)

    def test_decode_sequence_ids_length_exceeds_payload(self):
        header = SEQ_ID_MAGIC + bytes([SEQ_ID_VERSION])
        header += encode_varint(5) + b"\x00"
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(header)

    def test_decode_sequence_ids_missing_mode_byte(self):
        header = SEQ_ID_MAGIC + bytes([SEQ_ID_VERSION])
        header += encode_varint(0)
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(header)

    def test_decode_sequence_ids_unknown_compression_mode(self):
        header = SEQ_ID_MAGIC + bytes([SEQ_ID_VERSION])
        header += encode_varint(1) + bytes([9]) + b"x"
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(header)

    def test_decode_sequence_ids_unsupported_version(self):
        header = SEQ_ID_MAGIC + bytes([99]) + encode_varint(0)
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(header)

    def test_decode_sequence_ids_requires_zstd(self, monkeypatch):
        payload = b"data"
        buffer = self.make_sequence_buffer(1, payload)
        monkeypatch.setattr("clipkit.ecomp.reader._ZSTD_AVAILABLE", False)
        monkeypatch.setattr("clipkit.ecomp.reader.zstd", None)
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(buffer)

    def test_decode_sequence_ids_entry_exceeds_length(self):
        inner = encode_varint(1) + encode_varint(5) + b"abc"
        compressed = zlib.compress(inner)
        buffer = self.make_sequence_buffer(2, compressed)
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(buffer)

    def test_decode_sequence_ids_entry_exceeds_declared_block(self):
        payload = encode_varint(1) + encode_varint(5) + b"a"
        buffer = self.make_v1_buffer(payload)
        with pytest.raises(EcompDecodeError):
            _decode_sequence_ids(buffer)

    def test_decode_sequence_ids_zstd_mode(self, monkeypatch):
        inner = encode_varint(1) + encode_varint(3) + b"abc"

        class Dummy:
            def decompress(self, data):
                assert data == b"payload"
                return inner

        dummy_module = types.SimpleNamespace(ZstdDecompressor=lambda: Dummy())
        monkeypatch.setattr("clipkit.ecomp.reader._ZSTD_AVAILABLE", True)
        monkeypatch.setattr("clipkit.ecomp.reader.zstd", dummy_module)

        buffer = self.make_sequence_buffer(1, b"payload")
        ids, remaining = _decode_sequence_ids(buffer)
        assert ids == ["abc"]
        assert remaining == b""


class TestPermutationDecoding:
    def test_extract_permutation_chunk_roundtrip(self):
        permutation = [2, 0, 1]
        payload = bytearray(PERM_MAGIC)
        payload.append(1)  # version
        payload.append(0)  # flags: width=1, uncompressed
        payload.extend(encode_varint(len(permutation)))
        payload.extend(encode_varint(len(permutation)))
        payload.extend(bytes(permutation))
        raw = bytes(payload) + b"after"
        remaining, decoded = _extract_permutation_chunk(
            raw, {"length": len(payload)}
        )
        assert decoded == permutation
        assert remaining == b"after"

    def test_decode_permutation_from_list(self):
        assert _decode_permutation([3, 1, 2]) == [3, 1, 2]

    def test_decode_permutation_dict_with_zlib(self):
        permutation = [4, 0, 2]
        payload = bytes().join(
            int(value).to_bytes(1, "little") for value in permutation
        )
        compressed = zlib.compress(payload)
        metadata = {
            "version": 1,
            "size": len(permutation),
            "dtype": "uint8",
            "compression": "zlib",
            "data": base64.b64encode(compressed).decode("ascii"),
        }
        assert _decode_permutation(metadata) == permutation

    def test_decode_permutation_rejects_unknown_dtype(self):
        metadata = {
            "version": 1,
            "size": 1,
            "dtype": "uint64",
            "compression": "none",
            "data": "",
        }
        with pytest.raises(EcompDecodeError):
            _decode_permutation(metadata)

    def test_decode_permutation_rejects_payload_encoding(self):
        with pytest.raises(EcompDecodeError):
            _decode_permutation({"encoding": "payload"})

    def test_decode_permutation_rejects_format(self):
        with pytest.raises(EcompDecodeError):
            _decode_permutation("abc")

    def test_decode_permutation_bad_version(self):
        metadata = {
            "version": 3,
            "size": 1,
            "dtype": "uint8",
            "compression": "none",
            "data": "",
        }
        with pytest.raises(EcompDecodeError):
            _decode_permutation(metadata)

    def test_decode_permutation_bad_compression(self):
        metadata = {
            "version": 1,
            "size": 1,
            "dtype": "uint8",
            "compression": "bz2",
            "data": "",
        }
        with pytest.raises(EcompDecodeError):
            _decode_permutation(metadata)

    def test_decode_permutation_none_returns_empty(self):
        assert _decode_permutation(None) == []

    def test_decode_permutation_payload_length_mismatch(self):
        payload = base64.b64encode(b"ab").decode("ascii")
        metadata = {
            "version": 1,
            "size": 3,
            "dtype": "uint8",
            "compression": "none",
            "data": payload,
        }
        with pytest.raises(EcompDecodeError):
            _decode_permutation(metadata)


class TestPermutationChunk:
    def test_extract_permutation_invalid_length(self):
        with pytest.raises(EcompDecodeError):
            _extract_permutation_chunk(b"", {"length": 0})

    def test_extract_permutation_payload_short(self):
        with pytest.raises(EcompDecodeError):
            _extract_permutation_chunk(b"xx", {"length": 4})

    def test_extract_permutation_missing_magic(self):
        payload = b"BADA" + b"rest"
        with pytest.raises(EcompDecodeError):
            _extract_permutation_chunk(payload, {"length": len(payload)})

    def test_extract_permutation_bad_version(self):
        chunk = bytearray(PERM_MAGIC)
        chunk.append(9)
        chunk.append(0)
        chunk.extend(encode_varint(0))
        chunk.extend(encode_varint(0))
        with pytest.raises(EcompDecodeError):
            _extract_permutation_chunk(bytes(chunk), {"length": len(chunk)})

    def test_extract_permutation_bad_width_code(self):
        chunk = bytearray(PERM_MAGIC)
        chunk.append(PERM_VERSION)
        chunk.append(0x06)
        chunk.extend(encode_varint(0))
        chunk.extend(encode_varint(0))
        with pytest.raises(EcompDecodeError):
            _extract_permutation_chunk(bytes(chunk), {"length": len(chunk)})

    def test_extract_permutation_length_mismatch(self):
        payload = bytearray(PERM_MAGIC)
        payload.append(PERM_VERSION)
        payload.append(0)
        payload.extend(encode_varint(1))
        payload.extend(encode_varint(2))
        payload.extend(b"x")
        with pytest.raises(EcompDecodeError):
            _extract_permutation_chunk(bytes(payload), {"length": len(payload)})

    def test_extract_permutation_payload_size_mismatch(self):
        payload = bytearray(PERM_MAGIC)
        payload.append(PERM_VERSION)
        payload.append(0)
        payload.extend(encode_varint(3))
        payload.extend(encode_varint(2))
        payload.extend(b"ab")
        with pytest.raises(EcompDecodeError):
            _extract_permutation_chunk(bytes(payload), {"length": len(payload)})
class TestIterDeviationIndices:
    def test_iter_deviation_indices_matches_bitmask(self):
        bitmask = bytes([0b00001101])
        indices = _iter_deviation_indices(bitmask, 6)
        assert indices == [0, 2, 3]


class TestPackingHelpers:
    def test_pack_codes_pads_final_byte(self):
        result = _pack_codes([1, 3], 3)
        assert result == bytes([0b00101100])

    def test_unpack_codes_roundtrip(self):
        codes = [6, 1, 3, 4]
        data = _pack_codes(codes, 3)
        values = _unpack_codes(data, len(codes), 3)
        assert values == codes

    def test_read_varint_roundtrip(self):
        buffer = bytes([0x81, 0x01, 0x00])
        value, cursor = _read_varint(buffer, 0)
        assert (value, cursor) == (129, 2)
        value, cursor = _read_varint(buffer, cursor)
        assert (value, cursor) == (0, 3)

    def test_popcount_counts_bits(self):
        assert _popcount(bytes([0b11110000, 0b00001111])) == 8

    def test_read_varint_truncated(self):
        with pytest.raises(EcompDecodeError):
            _read_varint(bytes(), 0)

    def test_read_varint_too_long(self):
        data = bytes([0x80] * 9)
        with pytest.raises(EcompDecodeError):
            _read_varint(data, 0)

    def test_unpack_codes_zero_count(self):
        assert _unpack_codes(b"anything", 0, 1) == []

    def test_unpack_codes_truncated(self):
        with pytest.raises(EcompDecodeError):
            _unpack_codes(b"", 1, 1)

    def test_pack_codes_empty(self):
        assert _pack_codes([], 1) == b""

    def test_canonical_code_maps_skips_zero_lengths(self):
        encode_map, decode_map, max_len = _canonical_code_maps(["A", "B"], [0, 1])
        assert "A" not in encode_map
        assert max_len == 1


class TestDecodeVarint:
    def test_decode_varint_truncated(self):
        with pytest.raises(EcompDecodeError):
            _decode_varint(memoryview(b""), 0)

    def test_decode_varint_too_long(self):
        data = memoryview(bytes([0x80] * 9))
        with pytest.raises(EcompDecodeError):
            _decode_varint(data, 0)


class TestDecodeBlocks:
    def test_decode_blocks_empty_payload(self):
        assert decode_blocks(b"", 1, 1, []) == []

    def test_decode_blocks_dictionary_and_literal(self):
        alphabet = ["A", "C", "G"]
        payload = bytearray()
        payload.append(1)  # table count
        payload.append(ord("A"))
        payload.append(0)  # mode 0 table
        payload.append(2)  # entry count
        payload.extend(b"CG")
        payload.append(1)  # width

        dict_residues = _pack_codes([1], 1)
        payload.append(1)  # dictionary count
        payload.append(ord("A"))
        payload.append(0)  # copy mode
        payload.extend(encode_varint(1))  # deviation count
        mask_payload = bytes([0b00000001])
        payload.extend(encode_varint(len(mask_payload)))
        payload.extend(mask_payload)
        payload.extend(struct.pack(">H", len(dict_residues)))
        payload.extend(dict_residues)

        payload.extend(struct.pack(">I", 2))  # block count
        payload.append(1)  # dictionary marker
        payload.append(0)  # dictionary index
        payload.append(2)  # run length

        payload.append(0)  # literal marker
        payload.append(1)  # run length
        payload.append(ord("A"))
        payload.append(1)  # bitmask mode 1
        payload.extend(encode_varint(2))  # deviation count
        positions = encode_varint(1) + encode_varint(1)
        mask_payload = encode_varint(len(positions)) + positions
        payload.extend(encode_varint(len(mask_payload)))
        payload.extend(mask_payload)
        literal_residues = _pack_codes([0, 1], 1)
        payload.extend(struct.pack(">H", len(literal_residues)))
        payload.extend(literal_residues)

        blocks = decode_blocks(bytes(payload), 1, 2, alphabet)

        assert len(blocks) == 2
        assert blocks[0].run_length == 2
        assert blocks[1].run_length == 1
        assert blocks[1].bitmask != bytes([0])

    def test_decode_blocks_huffman_and_rle(self):
        alphabet = ["A", "C", "G", "T"]
        payload = bytearray()
        payload.append(1)  # table count
        payload.append(ord("C"))
        payload.append(1)  # huffman mode
        payload.append(2)  # entry count
        payload.extend(b"AT")
        payload.extend(bytes([1, 2]))

        payload.append(0)  # dictionary count
        payload.extend(struct.pack(">I", 1))

        payload.append(0)  # literal marker
        payload.append(1)  # run length
        payload.append(ord("C"))
        payload.append(2)  # bitmask mode 2
        payload.extend(encode_varint(1))  # deviation count
        mask_payload = bytes([0b00001000, 1])
        payload.extend(encode_varint(len(mask_payload)))
        payload.extend(mask_payload)

        encode_map, _, _ = _canonical_code_maps(["A", "T"], [1, 2])
        code, length = encode_map["T"]
        bits = f"{code:0{length}b}".ljust(8, "0")
        residue_bytes = int(bits, 2).to_bytes(1, "big")
        payload.extend(struct.pack(">H", len(residue_bytes)))
        payload.extend(residue_bytes)

        blocks = decode_blocks(bytes(payload), 1, 2, alphabet)

        assert len(blocks) == 1
        assert blocks[0].bitmask == bytes([0b00001000])

    def test_decode_blocks_consensus_table_truncated(self):
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes([1]), 1, 1, [])
        assert "Consensus table truncated" in str(exc.value)

    def test_decode_blocks_residues_truncated(self):
        payload = bytes([1, ord("A"), 0, 1])
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(payload, 1, 1, [])
        assert "Consensus residues truncated" in str(exc.value)

    def test_decode_blocks_width_missing(self):
        payload = bytes([1, ord("A"), 0, 1, ord("C")])
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(payload, 1, 1, [])
        assert "Consensus width missing" in str(exc.value)

    def test_decode_blocks_code_lengths_truncated(self):
        payload = bytes([1, ord("A"), 1, 1, ord("C")])
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(payload, 1, 1, [])
        assert "Consensus code lengths truncated" in str(exc.value)

    def test_decode_blocks_unknown_residue_mode(self):
        payload = bytes([1, ord("A"), 3, 0])
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(payload, 1, 1, [])
        assert "Unknown residue encoding mode" in str(exc.value)

    def test_decode_blocks_missing_dictionary_section(self):
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes([0]), 1, 1, [])
        assert "Missing dictionary section" in str(exc.value)

    def test_decode_blocks_dictionary_entry_truncated(self):
        payload = bytes([0, 1])
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(payload, 1, 1, [])
        assert "Dictionary entry truncated" in str(exc.value)

    def test_decode_blocks_dictionary_mask_truncated(self):
        payload = bytearray([0, 1, ord("A"), 0])
        payload.extend(encode_varint(1))
        payload.extend(encode_varint(2))
        payload.append(1)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Dictionary mask truncated" in str(exc.value)

    def test_decode_blocks_dictionary_residue_header_truncated(self):
        payload = bytearray([0, 1, ord("A"), 0])
        payload.extend(encode_varint(1))
        payload.extend(encode_varint(1))
        payload.append(1)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Dictionary residue header truncated" in str(exc.value)

    def test_decode_blocks_dictionary_residues_truncated(self):
        payload = bytearray([0, 1, ord("A"), 0])
        payload.extend(encode_varint(1))
        payload.extend(encode_varint(0))
        payload.extend(struct.pack(">H", 2))
        payload.append(0)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Dictionary residues truncated" in str(exc.value)

    def test_decode_blocks_missing_block_count(self):
        payload = bytearray([0, 0])
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Missing block count" in str(exc.value)

    def test_decode_blocks_block_data_truncated(self):
        payload = bytearray([0, 0])
        payload.extend(struct.pack(">I", 1))
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Block data truncated" in str(exc.value)

    def test_decode_blocks_dictionary_block_truncated(self):
        payload = bytearray([0, 0])
        payload.extend(struct.pack(">I", 1))
        payload.append(1)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Dictionary block truncated" in str(exc.value)

    def test_decode_blocks_dictionary_index_out_of_range(self):
        payload = bytearray([0, 1, ord("A"), 0])
        payload.extend(encode_varint(0))
        payload.extend(encode_varint(0))
        payload.extend(struct.pack(">H", 0))
        payload.extend(struct.pack(">I", 1))
        payload.append(1)
        payload.append(3)
        payload.append(1)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Dictionary index out of range" in str(exc.value)

    def test_decode_blocks_literal_block_truncated(self):
        payload = bytearray([0, 0])
        payload.extend(struct.pack(">I", 1))
        payload.append(0)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Literal block truncated" in str(exc.value)

    def test_decode_blocks_literal_mask_truncated(self):
        payload = bytearray([0, 0])
        payload.extend(struct.pack(">I", 1))
        payload.extend([0, 1, ord("A"), 0])
        payload.extend(encode_varint(1))
        payload.extend(encode_varint(2))
        payload.append(1)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Literal mask truncated" in str(exc.value)

    def test_decode_blocks_literal_residue_header_truncated(self):
        payload = bytearray([0, 0])
        payload.extend(struct.pack(">I", 1))
        payload.extend([0, 1, ord("A"), 0])
        payload.extend(encode_varint(1))
        payload.extend(encode_varint(0))
        payload.append(0)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Literal residue header truncated" in str(exc.value)

    def test_decode_blocks_literal_residues_truncated(self):
        payload = bytearray([0, 0])
        payload.extend(struct.pack(">I", 1))
        payload.extend([0, 1, ord("A"), 0])
        payload.extend(encode_varint(1))
        payload.extend(encode_varint(0))
        payload.extend(struct.pack(">H", 2))
        payload.append(1)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Literal residues truncated" in str(exc.value)

    def test_decode_blocks_unknown_block_marker(self):
        payload = bytearray([0, 0])
        payload.extend(struct.pack(">I", 1))
        payload.append(5)
        with pytest.raises(EcompDecodeError) as exc:
            decode_blocks(bytes(payload), 1, 1, [])
        assert "Unknown block marker" in str(exc.value)

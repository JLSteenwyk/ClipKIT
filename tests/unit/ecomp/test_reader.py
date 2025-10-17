import base64
import json
import zlib

import pytest

from clipkit.ecomp.reader import (
    METADATA_COMPRESSED_MAGIC,
    METADATA_CODEC_VERSION,
    PERM_MAGIC,
    SEQ_ID_MAGIC,
    SEQ_ID_VERSION,
    EcompDecodeError,
    _canonical_code_maps,
    _decode_bitmask,
    _decode_metadata,
    _decode_permutation,
    _decode_residue_stream,
    _decode_residues,
    _decode_sequence_ids,
    _extract_permutation_chunk,
    _iter_deviation_indices,
    _pack_codes,
    _popcount,
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

    def test_decode_residues_roundtrip(self):
        data = _pack_codes([0, 2, 1], 2)
        alphabet = ["A", "C", "G"]
        result = _decode_residues(data, 3, 2, alphabet)
        assert result == ["A", "G", "C"]

    def test_decode_residues_insufficient_data(self):
        data = bytes([0b10000000])
        with pytest.raises(EcompDecodeError):
            _decode_residues(data, 2, 6, ["A"])


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

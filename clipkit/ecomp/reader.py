"""Read `.ecomp` archives into Biopython alignments."""

from __future__ import annotations

import base64
import gzip
import io
import json
import lzma
import math
import struct
import zlib
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Iterable, Sequence

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:  # pragma: no cover - optional dependency
    import zstandard as zstd

    _ZSTD_AVAILABLE = True
except ModuleNotFoundError:  # pragma: no cover - zstandard not installed
    zstd = None
    _ZSTD_AVAILABLE = False

HEADER_MAGIC = b"ECOMP001"
LEGACY_HEADER_STRUCT = ">8sBBBQ"
HEADER_STRUCT = ">8sBBBQQ"
INLINE_METADATA_VERSION = (0, 2, 0)
METADATA_SUFFIX = ".json"
METADATA_COMPRESSED_MAGIC = b"ECMZ"
METADATA_CODEC_VERSION = 1

SEQ_ID_MAGIC = b"ECID"
SEQ_ID_VERSION = 2
PERM_MAGIC = b"ECPE"
PERM_VERSION = 1


class EcompDecodeError(RuntimeError):
    """Raised when `.ecomp` parsing encounters invalid data."""


@dataclass
class RunLengthBlock:
    """A run of identical consensus columns encoded with bitmasks."""

    consensus: str
    bitmask: bytes
    residues: bytes
    run_length: int


def read_ecomp(
    path: str | Path,
    *,
    metadata_path: str | Path | None = None,
) -> tuple[MultipleSeqAlignment, dict[str, object]]:
    """Return ``(alignment, metadata)`` decoded from an `.ecomp` archive."""

    payload, metadata, version = _read_archive(
        path, metadata_path=metadata_path
    )
    alignment = _decompress_alignment(payload, metadata)
    alignment.annotations.setdefault("ecomp_metadata", metadata)
    alignment.annotations.setdefault("ecomp_version", version)
    return alignment, metadata


def _read_archive(
    path: str | Path,
    *,
    metadata_path: str | Path | None = None,
) -> tuple[bytes, dict[str, object], tuple[int, int, int]]:
    file_path = Path(path)
    data = file_path.read_bytes()
    legacy_size = struct.calcsize(LEGACY_HEADER_STRUCT)
    if len(data) < legacy_size:
        raise EcompDecodeError("File too small to be a valid .ecomp archive")

    header = data[:legacy_size]
    magic, major, minor, patch, payload_len = struct.unpack(
        LEGACY_HEADER_STRUCT, header
    )
    if magic != HEADER_MAGIC:
        raise EcompDecodeError("Invalid .ecomp magic header")

    version = (major, minor, patch)
    cursor = legacy_size

    metadata_len: int | None = None
    if version >= INLINE_METADATA_VERSION:
        extended_size = struct.calcsize(HEADER_STRUCT)
        if len(data) < extended_size:
            raise EcompDecodeError("Archive truncated before metadata header")
        (metadata_len,) = struct.unpack(
            ">Q", data[legacy_size:extended_size]
        )
        cursor = extended_size

    payload_end = cursor + payload_len
    if len(data) < payload_end:
        raise EcompDecodeError("Payload length exceeds available data")
    payload = data[cursor:payload_end]

    if metadata_len is not None:
        metadata_end = payload_end + metadata_len
        if len(data) < metadata_end:
            raise EcompDecodeError("Metadata length exceeds available data")
        metadata_bytes = data[payload_end:metadata_end]
        metadata = _decode_metadata(metadata_bytes)
    else:
        meta_path = (
            Path(metadata_path)
            if metadata_path
            else _derive_metadata_path(file_path)
        )
        if not meta_path.exists():
            raise FileNotFoundError(
                "Legacy .ecomp archives require metadata sidecar"
            )
        metadata = _read_metadata(meta_path)

    return payload, metadata, version


def _decode_metadata(data: bytes) -> dict[str, object]:
    if data.startswith(METADATA_COMPRESSED_MAGIC):
        header_len = len(METADATA_COMPRESSED_MAGIC)
        if len(data) <= header_len:
            raise EcompDecodeError("Compressed metadata header truncated")
        codec_version = data[header_len]
        if codec_version != METADATA_CODEC_VERSION:
            raise EcompDecodeError(
                f"Unsupported compressed metadata version: {codec_version}"
            )
        payload = zlib.decompress(data[header_len + 1 :])
    else:
        payload = data
    try:
        return json.loads(payload.decode("utf-8"))
    except json.JSONDecodeError as exc:
        raise EcompDecodeError("Metadata JSON is malformed") from exc


def _read_metadata(path: Path) -> dict[str, object]:
    return _decode_metadata(path.read_bytes())


def _derive_metadata_path(ecomp_path: Path) -> Path:
    return ecomp_path.with_suffix(METADATA_SUFFIX)


def _decompress_alignment(
    payload: bytes,
    metadata: dict[str, object],
) -> MultipleSeqAlignment:
    fallback = metadata.get("fallback")
    if isinstance(fallback, dict) and fallback.get("type") == "gzip":
        fasta_bytes = gzip.decompress(payload)
        alignment = _alignment_from_fasta_bytes(
            fasta_bytes, metadata.get("source_format", "fasta")
        )
        return alignment

    alignment_length = int(metadata["alignment_length"])
    num_sequences = int(metadata["num_sequences"])

    sequence_ids = metadata.get("sequence_ids")
    if sequence_ids is not None and len(sequence_ids) != num_sequences:
        raise EcompDecodeError(
            "Sequence ID count does not match `num_sequences` in metadata"
        )

    alphabet = metadata.get("alphabet", [])
    bitmask_bytes = metadata.get("bitmask_bytes")
    bits_per_symbol = metadata.get("bits_per_symbol")

    if not isinstance(bitmask_bytes, int) or bitmask_bytes <= 0:
        raise EcompDecodeError("Metadata missing valid `bitmask_bytes`")
    if not isinstance(bits_per_symbol, int) or bits_per_symbol <= 0:
        if alphabet:
            bits_per_symbol = max(1, math.ceil(math.log2(len(alphabet))))
        else:
            raise EcompDecodeError("Metadata missing valid `bits_per_symbol`")

    encoding = metadata.get("payload_encoding", "raw")
    if encoding == "zstd":
        if not _ZSTD_AVAILABLE:
            raise EcompDecodeError(
                "Payload encoded with zstd but `zstandard` is not installed"
            )
        payload_bytes = zstd.ZstdDecompressor().decompress(payload)
    elif encoding == "xz":
        payload_bytes = lzma.decompress(payload)
    elif encoding == "zlib":
        payload_bytes = zlib.decompress(payload)
    elif encoding in {"raw", None}:
        payload_bytes = payload
    else:
        raise EcompDecodeError(f"Unsupported payload encoding: {encoding}")

    permutation: list[int] | None = None
    perm_meta = metadata.get("sequence_permutation")
    if isinstance(perm_meta, dict) and perm_meta.get("encoding") == "payload":
        payload_bytes, permutation = _extract_permutation_chunk(
            payload_bytes, perm_meta
        )
        metadata["sequence_permutation"] = permutation
    elif perm_meta is not None:
        permutation = _decode_permutation(perm_meta)
        metadata["sequence_permutation"] = permutation

    if payload_bytes.startswith(SEQ_ID_MAGIC):
        seq_ids, payload_bytes = _decode_sequence_ids(payload_bytes)
        if sequence_ids is None:
            sequence_ids = seq_ids
            metadata["sequence_ids"] = seq_ids
        elif list(sequence_ids) != seq_ids:
            raise EcompDecodeError(
                "Sequence IDs in payload do not match metadata"
            )
    if sequence_ids is None:
        raise EcompDecodeError("Sequence identifiers missing from metadata")

    blocks = decode_blocks(
        payload_bytes, bitmask_bytes, bits_per_symbol, alphabet
    )
    sequences = [
        ["" for _ in range(alignment_length)] for _ in range(num_sequences)
    ]
    alphabet_list = list(alphabet)

    column = 0
    for block in blocks:
        for _ in range(block.run_length):
            if column >= alignment_length:
                raise EcompDecodeError(
                    "Decoded columns exceed expected alignment length"
                )
            indices = _iter_deviation_indices(
                block.bitmask, num_sequences
            )
            residues = _decode_residues(
                block.residues,
                len(indices),
                bits_per_symbol,
                alphabet_list,
            )
            for row in sequences:
                row[column] = block.consensus
            for seq_idx, residue in zip(indices, residues):
                sequences[seq_idx][column] = residue
            column += 1

    if column != alignment_length:
        raise EcompDecodeError(
            "Decoded columns do not match expected alignment length"
        )

    reconstructed = ["".join(row) for row in sequences]

    if permutation:
        inverse = [0] * len(permutation)
        for new_pos, original_index in enumerate(permutation):
            inverse[original_index] = new_pos
        reordered = [
            reconstructed[inverse[idx]] for idx in range(len(inverse))
        ]
        id_reordered = [
            sequence_ids[inverse[idx]] for idx in range(len(inverse))
        ]
        reconstructed = reordered
        sequence_ids = id_reordered

    _validate_checksum(reconstructed, metadata.get("checksum_sha256"))

    records = []
    for seq_id, seq in zip(sequence_ids, reconstructed):
        seq_id_str = str(seq_id)
        records.append(
            SeqRecord(Seq(seq), id=seq_id_str, description=seq_id_str)
        )
    alignment = MultipleSeqAlignment(records)
    alignment.annotations["source_format"] = metadata.get(
        "source_format", "ecomp"
    )
    return alignment


def _validate_checksum(sequences: Iterable[str], checksum: object) -> None:
    if not checksum:
        return
    digest = sha256()
    for seq in sequences:
        digest.update(seq.encode("utf-8"))
        digest.update(b"\n")
    if digest.hexdigest() != checksum:
        raise EcompDecodeError(
            "Checksum mismatch while decoding .ecomp archive"
        )


def _alignment_from_fasta_bytes(
    data: bytes, source_format: str | None
) -> MultipleSeqAlignment:
    buffer = io.StringIO(data.decode("utf-8"))
    records: list[SeqRecord] = []
    current_id: str | None = None
    current_seq: list[str] = []
    for line in buffer:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if current_id is not None:
                seq_value = "".join(current_seq)
                records.append(
                    SeqRecord(
                        Seq(seq_value),
                        id=current_id,
                        description=current_id,
                    )
                )
            current_id = stripped[1:]
            current_seq = []
        else:
            current_seq.append(stripped)
    if current_id is not None:
        seq_value = "".join(current_seq)
        records.append(
            SeqRecord(
                Seq(seq_value), id=current_id, description=current_id
            )
        )
    alignment = MultipleSeqAlignment(records)
    alignment.annotations["source_format"] = source_format or "fasta"
    return alignment


def decode_blocks(
    payload: bytes,
    bitmask_bytes: int,
    bits_per_symbol: int,
    alphabet: Sequence[str],
) -> list[RunLengthBlock]:
    blocks: list[RunLengthBlock] = []
    cursor = 0
    total = len(payload)
    if cursor >= total:
        return blocks

    table_count = payload[cursor]
    cursor += 1
    models: dict[str, dict[str, object]] = {}
    for _ in range(table_count):
        if cursor + 3 > total:
            raise EcompDecodeError("Consensus table truncated")
        consensus_value = payload[cursor]
        cursor += 1
        mode = payload[cursor]
        cursor += 1
        entry_count = payload[cursor]
        cursor += 1
        residues: list[str] = []
        for _ in range(entry_count):
            if cursor >= total:
                raise EcompDecodeError("Consensus residues truncated")
            residues.append(bytes([payload[cursor]]).decode("ascii"))
            cursor += 1
        if mode == 0:
            if cursor >= total:
                raise EcompDecodeError("Consensus width missing")
            width = payload[cursor]
            cursor += 1
            models[chr(consensus_value)] = {
                "mode": 0,
                "residues": residues,
                "width": max(1, width),
                "decode_map": {idx: res for idx, res in enumerate(residues)},
            }
        elif mode == 1:
            lengths: list[int] = []
            for _ in range(entry_count):
                if cursor >= total:
                    raise EcompDecodeError("Consensus code lengths truncated")
                lengths.append(payload[cursor])
                cursor += 1
            encode_map, decode_map, max_len = _canonical_code_maps(
                residues, lengths
            )
            models[chr(consensus_value)] = {
                "mode": 1,
                "residues": residues,
                "lengths": lengths,
                "decode_map": decode_map,
                "max_code_len": max_len,
            }
        else:
            raise EcompDecodeError("Unknown residue encoding mode")

    if cursor >= total:
        raise EcompDecodeError("Missing dictionary section")

    dict_count = payload[cursor]
    cursor += 1
    dictionary: list[tuple[str, bytes, bytes]] = []
    for _ in range(dict_count):
        if cursor + 2 > total:
            raise EcompDecodeError("Dictionary entry truncated")
        consensus_value = payload[cursor]
        cursor += 1
        mode = payload[cursor]
        cursor += 1
        deviation_count, cursor = _read_varint(payload, cursor)
        mask_len, cursor = _read_varint(payload, cursor)
        if cursor + mask_len > total:
            raise EcompDecodeError("Dictionary mask truncated")
        mask_payload = payload[cursor : cursor + mask_len]
        cursor += mask_len
        if cursor + 2 > total:
            raise EcompDecodeError("Dictionary residue header truncated")
        residues_len = struct.unpack(">H", payload[cursor : cursor + 2])[0]
        cursor += 2
        if cursor + residues_len > total:
            raise EcompDecodeError("Dictionary residues truncated")
        residues = payload[cursor : cursor + residues_len]
        cursor += residues_len
        bitmask = _decode_bitmask(
            mode, mask_payload, deviation_count, bitmask_bytes
        )
        dictionary.append((chr(consensus_value), bitmask, residues))

    if cursor + 4 > total:
        raise EcompDecodeError("Missing block count")
    (block_count,) = struct.unpack(">I", payload[cursor : cursor + 4])
    cursor += 4

    alphabet_lookup = {char: idx for idx, char in enumerate(alphabet)}

    for _ in range(block_count):
        if cursor >= total:
            raise EcompDecodeError("Block data truncated")
        marker = payload[cursor]
        cursor += 1
        if marker == 1:
            if cursor + 2 > total:
                raise EcompDecodeError("Dictionary block truncated")
            dict_id = payload[cursor]
            cursor += 1
            run_length = payload[cursor]
            cursor += 1
            try:
                consensus, bitmask, residues = dictionary[dict_id]
            except IndexError as exc:
                raise EcompDecodeError(
                    "Dictionary index out of range"
                ) from exc
            decoded = _decode_residue_stream(
                consensus,
                bitmask,
                residues,
                models,
                bits_per_symbol,
                alphabet_lookup,
            )
            blocks.append(
                RunLengthBlock(
                    consensus=consensus,
                    bitmask=bitmask,
                    residues=decoded,
                    run_length=run_length,
                )
            )
        elif marker == 0:
            if cursor + 3 > total:
                raise EcompDecodeError("Literal block truncated")
            run_length = payload[cursor]
            cursor += 1
            consensus_value = payload[cursor]
            cursor += 1
            mode = payload[cursor]
            cursor += 1
            deviation_count, cursor = _read_varint(payload, cursor)
            mask_len, cursor = _read_varint(payload, cursor)
            if cursor + mask_len > total:
                raise EcompDecodeError("Literal mask truncated")
            mask_payload = payload[cursor : cursor + mask_len]
            cursor += mask_len
            if cursor + 2 > total:
                raise EcompDecodeError("Literal residue header truncated")
            residues_len = struct.unpack(">H", payload[cursor : cursor + 2])[0]
            cursor += 2
            if cursor + residues_len > total:
                raise EcompDecodeError("Literal residues truncated")
            residues = payload[cursor : cursor + residues_len]
            cursor += residues_len
            bitmask = _decode_bitmask(
                mode, mask_payload, deviation_count, bitmask_bytes
            )
            consensus = chr(consensus_value)
            decoded = _decode_residue_stream(
                consensus,
                bitmask,
                residues,
                models,
                bits_per_symbol,
                alphabet_lookup,
            )
            blocks.append(
                RunLengthBlock(
                    consensus=consensus,
                    bitmask=bitmask,
                    residues=decoded,
                    run_length=run_length,
                )
            )
        else:
            raise EcompDecodeError("Unknown block marker")
    return blocks


def _read_varint(data: bytes, cursor: int) -> tuple[int, int]:
    shift = 0
    result = 0
    while True:
        if cursor >= len(data):
            raise EcompDecodeError(
                "Unexpected end of data while reading varint"
            )
        byte = data[cursor]
        cursor += 1
        result |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            return result, cursor
        shift += 7
        if shift > 56:
            raise EcompDecodeError("Varint too long")


def _decode_bitmask(
    mode: int,
    payload: bytes,
    deviation_count: int,
    bitmask_bytes: int,
) -> bytes:
    if deviation_count == 0:
        return bytes(bitmask_bytes)
    if mode == 0:
        bitmask = bytearray(bitmask_bytes)
        length = min(len(payload), bitmask_bytes)
        bitmask[:length] = payload[:length]
        return bytes(bitmask)
    if mode == 1:
        length, cursor = _read_varint(payload, 0)
        encoded = payload[cursor : cursor + length]
        positions = _decode_positions(encoded, deviation_count)
        bitmask = bytearray(bitmask_bytes)
        max_bits = bitmask_bytes * 8
        for pos in positions:
            if pos >= max_bits:
                raise EcompDecodeError(
                    "Bitmask position exceeds sequence count"
                )
            byte_index = pos // 8
            bit_index = pos % 8
            bitmask[byte_index] |= 1 << bit_index
        return bytes(bitmask)
    if mode == 2:
        bitmask = bytearray()
        idx = 0
        while idx < len(payload) and len(bitmask) < bitmask_bytes:
            if idx + 2 > len(payload):
                raise EcompDecodeError("RLE bitmask truncated")
            byte = payload[idx]
            count = payload[idx + 1]
            idx += 2
            bitmask.extend([byte] * count)
        if len(bitmask) < bitmask_bytes:
            bitmask.extend([0] * (bitmask_bytes - len(bitmask)))
        return bytes(bitmask[:bitmask_bytes])
    raise EcompDecodeError(f"Unknown bitmask encoding mode: {mode}")


def _decode_positions(encoded: bytes, count: int) -> list[int]:
    positions: list[int] = []
    cursor = 0
    previous = -1
    for _ in range(count):
        value, cursor = _read_varint(encoded, cursor)
        pos = previous + value
        positions.append(pos)
        previous = pos
    return positions


def _decode_residue_stream(
    consensus: str,
    bitmask: bytes,
    encoded: bytes,
    models: dict[str, dict[str, object]],
    bits_per_symbol: int,
    alphabet_lookup: dict[str, int],
) -> bytes:
    deviation_count = _popcount(bitmask)
    if deviation_count == 0:
        return b""
    model = models.get(consensus)
    if model is None:
        raise EcompDecodeError("Missing residue model for consensus")
    if model["mode"] == 0:
        width = model["width"]
        residues = model["residues"]
        local_codes = _unpack_codes(encoded, deviation_count, width)
        try:
            chars = [residues[code] for code in local_codes]
        except IndexError as exc:
            raise EcompDecodeError("Residue code exceeds table size") from exc
    else:
        decode_map = model["decode_map"]
        max_len = model["max_code_len"]
        chars: list[str] = []
        total_bits = len(encoded) * 8
        bit_index = 0
        for _ in range(deviation_count):
            current = 0
            length = 0
            while True:
                if bit_index >= total_bits:
                    raise EcompDecodeError(
                        "Insufficient Huffman bits for residue stream"
                    )
                byte = encoded[bit_index // 8]
                bit = (byte >> (7 - (bit_index % 8))) & 1
                bit_index += 1
                current = (current << 1) | bit
                length += 1
                residue = decode_map.get((length, current))
                if residue is not None:
                    chars.append(residue)
                    break
                if length > max_len:
                    raise EcompDecodeError(
                        "Invalid Huffman code in residue stream"
                    )
    try:
        global_codes = [alphabet_lookup[char] for char in chars]
    except KeyError as exc:
        raise EcompDecodeError(
            "Residue not present in alphabet lookup"
        ) from exc
    return _pack_codes(global_codes, bits_per_symbol)


def _unpack_codes(
    data: bytes, count: int, bits_per_symbol: int
) -> list[int]:
    if count == 0:
        return []
    values: list[int] = []
    buffer = 0
    bits_in_buffer = 0
    mask = (1 << bits_per_symbol) - 1
    for byte in data:
        buffer = (buffer << 8) | byte
        bits_in_buffer += 8
        while bits_in_buffer >= bits_per_symbol and len(values) < count:
            bits_in_buffer -= bits_per_symbol
            values.append((buffer >> bits_in_buffer) & mask)
            buffer &= (1 << bits_in_buffer) - 1
        if len(values) == count:
            break
    if len(values) != count:
        raise EcompDecodeError("Residue stream truncated while unpacking")
    return values


def _pack_codes(codes: Sequence[int], bits_per_symbol: int) -> bytes:
    if not codes:
        return b""
    buffer = 0
    bits_in_buffer = 0
    output = bytearray()
    for code in codes:
        buffer = (buffer << bits_per_symbol) | code
        bits_in_buffer += bits_per_symbol
        while bits_in_buffer >= 8:
            bits_in_buffer -= 8
            output.append((buffer >> bits_in_buffer) & 0xFF)
            buffer &= (1 << bits_in_buffer) - 1
    if bits_in_buffer:
        output.append((buffer << (8 - bits_in_buffer)) & 0xFF)
    return bytes(output)


def _canonical_code_maps(
    residues: Sequence[str], lengths: Sequence[int]
) -> tuple[dict[str, tuple[int, int]], dict[tuple[int, int], str], int]:
    items = sorted(zip(residues, lengths), key=lambda item: (item[1], item[0]))
    encode_map: dict[str, tuple[int, int]] = {}
    decode_map: dict[tuple[int, int], str] = {}
    code = 0
    prev_len = 0
    max_len = 0
    for residue, length in items:
        if length <= 0:
            continue
        code <<= length - prev_len
        encode_map[residue] = (code, length)
        decode_map[(length, code)] = residue
        prev_len = length
        code += 1
        if length > max_len:
            max_len = length
    return encode_map, decode_map, max_len


def _popcount(data: bytes) -> int:
    return sum(bin(byte).count("1") for byte in data)


def _decode_sequence_ids(buffer: bytes) -> tuple[list[str], bytes]:
    view = memoryview(buffer)
    if len(view) < len(SEQ_ID_MAGIC) + 1:
        raise EcompDecodeError("Sequence ID block truncated")
    if not view.tobytes().startswith(SEQ_ID_MAGIC):
        raise EcompDecodeError("Sequence ID magic missing")
    cursor = len(SEQ_ID_MAGIC)
    version = view[cursor]
    cursor += 1
    block_length, cursor = _decode_varint(view, cursor)
    end = cursor + block_length
    if end > len(view):
        raise EcompDecodeError("Sequence ID block length exceeds payload")

    if version == 1:
        count, cursor = _decode_varint(view, cursor)
        ids: list[str] = []
        for _ in range(count):
            name_len, cursor = _decode_varint(view, cursor)
            if cursor + name_len > end:
                raise EcompDecodeError(
                    "Sequence ID entry exceeds declared block length"
                )
            name_bytes = view[cursor : cursor + name_len].tobytes()
            cursor += name_len
            ids.append(name_bytes.decode("utf-8"))
    elif version == SEQ_ID_VERSION:
        if cursor >= end:
            raise EcompDecodeError("Sequence ID block missing mode byte")
        mode = view[cursor]
        cursor += 1
        encoded = view[cursor:end].tobytes()
        if mode == 0:
            output = encoded
        elif mode == 1:
            if not _ZSTD_AVAILABLE:
                raise EcompDecodeError(
                    "Sequence IDs compressed with zstd but `zstandard` missing"
                )
            output = zstd.ZstdDecompressor().decompress(encoded)
        elif mode == 2:
            output = zlib.decompress(encoded)
        else:
            raise EcompDecodeError("Unsupported sequence ID compression mode")

        data_view = memoryview(output)
        cursor_out = 0
        count, cursor_out = _decode_varint(data_view, cursor_out)
        ids = []
        for _ in range(count):
            name_len, cursor_out = _decode_varint(data_view, cursor_out)
            if cursor_out + name_len > len(data_view):
                raise EcompDecodeError(
                    "Sequence ID entry exceeds decoded block length"
                )
            end_pos = cursor_out + name_len
            name_bytes = data_view[cursor_out:end_pos].tobytes()
            cursor_out += name_len
            ids.append(name_bytes.decode("utf-8"))
        cursor = end
    else:
        raise EcompDecodeError(
            f"Unsupported sequence ID block version: {version}"
        )

    if cursor != end:
        raise EcompDecodeError("Sequence ID block contains trailing data")

    remaining = view[end:].tobytes()
    return ids, remaining


def _decode_varint(buffer: memoryview, cursor: int) -> tuple[int, int]:
    shift = 0
    result = 0
    while True:
        if cursor >= len(buffer):
            raise EcompDecodeError(
                "Unexpected end of data while decoding varint"
            )
        byte = buffer[cursor]
        cursor += 1
        result |= (byte & 0x7F) << shift
        if not (byte & 0x80):
            return result, cursor
        shift += 7
        if shift > 56:
            raise EcompDecodeError("Varint too long")


def _extract_permutation_chunk(
    payload: bytes, metadata: dict[str, object]
) -> tuple[bytes, list[int]]:
    length = metadata.get("length")
    if not isinstance(length, int) or length <= 0:
        raise EcompDecodeError("Invalid permutation metadata length")
    if len(payload) < length:
        raise EcompDecodeError("Permutation chunk exceeds payload size")

    chunk = memoryview(payload[:length])
    cursor = 0
    if chunk[cursor : cursor + 4].tobytes() != PERM_MAGIC:
        raise EcompDecodeError("Permutation chunk missing magic header")
    cursor += 4

    version = chunk[cursor]
    cursor += 1
    if version != PERM_VERSION:
        raise EcompDecodeError(
            f"Unsupported permutation chunk version: {version}"
        )

    flags = chunk[cursor]
    cursor += 1
    compression_flag = flags & 0x01
    width_code = (flags >> 1) & 0x03
    width_lookup = {0: 1, 1: 2, 2: 4}
    try:
        width = width_lookup[width_code]
    except KeyError as exc:
        raise EcompDecodeError("Unsupported permutation width code") from exc

    size, cursor = _decode_varint(chunk, cursor)
    payload_len, cursor = _decode_varint(chunk, cursor)
    end = cursor + payload_len
    if end != length:
        raise EcompDecodeError("Permutation chunk length mismatch")

    payload_bytes = bytes(chunk[cursor:end])
    if compression_flag:
        payload_bytes = zlib.decompress(payload_bytes)

    if len(payload_bytes) != size * width:
        raise EcompDecodeError("Permutation payload size mismatch")

    permutation = [
        int.from_bytes(payload_bytes[i : i + width], "little", signed=False)
        for i in range(0, len(payload_bytes), width)
    ]

    remaining = payload[length:]
    return remaining, permutation


def _decode_permutation(value: object) -> list[int]:
    if value is None:
        return []
    if isinstance(value, list):
        return [int(v) for v in value]
    if not isinstance(value, dict):
        raise EcompDecodeError("Unsupported permutation metadata format")

    if value.get("encoding") == "payload":
        raise EcompDecodeError(
            "Payload-encoded permutation must be decoded from the payload"
        )

    version = value.get("version", 0)
    if version != PERM_VERSION:
        raise EcompDecodeError(
            f"Unsupported permutation metadata version: {version}"
        )

    size = int(value.get("size", 0))
    dtype = value.get("dtype")
    compression = value.get("compression", "none")
    data = value.get("data", "")

    payload = base64.b64decode(data.encode("ascii")) if data else b""
    if compression == "zlib" and payload:
        payload = zlib.decompress(payload)
    elif compression not in {"none", "zlib"}:
        raise EcompDecodeError(
            f"Unsupported permutation compression: {compression}"
        )

    width_lookup = {"uint8": 1, "uint16": 2, "uint32": 4}
    try:
        width = width_lookup[dtype]
    except KeyError as exc:
        raise EcompDecodeError("Unsupported permutation dtype") from exc

    expected = size * width
    if len(payload) != expected:
        raise EcompDecodeError("Permutation payload length mismatch")

    return [
        int.from_bytes(payload[i : i + width], "little", signed=False)
        for i in range(0, expected, width)
    ]


def _iter_deviation_indices(bitmask: bytes, num_sequences: int) -> list[int]:
    indices: list[int] = []
    for seq_idx in range(num_sequences):
        byte_index = seq_idx // 8
        bit_index = seq_idx % 8
        if bitmask[byte_index] & (1 << bit_index):
            indices.append(seq_idx)
    return indices


def _decode_residues(
    data: bytes,
    count: int,
    bits_per_symbol: int,
    alphabet: Sequence[str],
) -> list[str]:
    if count == 0:
        return []
    mask = (1 << bits_per_symbol) - 1
    values: list[int] = []
    buffer = 0
    bits_in_buffer = 0
    data_iter = iter(data)
    while len(values) < count:
        while bits_in_buffer < bits_per_symbol:
            try:
                byte = next(data_iter)
            except StopIteration as exc:
                raise EcompDecodeError(
                    "Insufficient residue data during decode"
                ) from exc
            buffer = (buffer << 8) | byte
            bits_in_buffer += 8
        shift = bits_in_buffer - bits_per_symbol
        value = (buffer >> shift) & mask
        buffer &= (1 << shift) - 1
        bits_in_buffer -= bits_per_symbol
        values.append(value)
    alphabet_list = list(alphabet)
    try:
        return [alphabet_list[value] for value in values]
    except IndexError as exc:
        raise EcompDecodeError("Residue code exceeds alphabet size") from exc

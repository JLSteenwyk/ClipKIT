"""Write Evolutionary Compression (.ecomp) archives."""

from __future__ import annotations

import gzip
import json
import struct
import zlib
from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Iterable

from Bio.Align import MultipleSeqAlignment

from .reader import (
    HEADER_MAGIC,
    HEADER_STRUCT,
    INLINE_METADATA_VERSION,
    METADATA_COMPRESSED_MAGIC,
    METADATA_CODEC_VERSION,
)


class EcompEncodeError(RuntimeError):
    """Raised when `.ecomp` encoding encounters invalid data."""


@dataclass
class EncodedArchive:
    payload: bytes
    metadata_bytes: bytes
    metadata: dict[str, object]
    version_tuple: tuple[int, int, int]


def write_ecomp(
    alignment: MultipleSeqAlignment,
    output_path: str | Path,
    base_metadata: dict[str, object] | None = None,
) -> dict[str, object]:
    """Serialize *alignment* to ``output_path`` as an `.ecomp` archive.

    The encoder always emits a gzip fallback payload so downstream tools can
    decode the result without the native ClipKIT structures.
    """

    archive = _build_archive(alignment, base_metadata)
    path = Path(output_path)
    header = struct.pack(
        HEADER_STRUCT,
        HEADER_MAGIC,
        archive.version_tuple[0],
        archive.version_tuple[1],
        archive.version_tuple[2],
        len(archive.payload),
        len(archive.metadata_bytes),
    )

    with path.open("wb") as handle:
        handle.write(header)
        handle.write(archive.payload)
        handle.write(archive.metadata_bytes)

    return archive.metadata


def _build_archive(
    alignment: MultipleSeqAlignment,
    base_metadata: dict[str, object] | None,
) -> EncodedArchive:
    fasta_bytes = _alignment_to_fasta_bytes(alignment)
    payload = gzip.compress(fasta_bytes)
    metadata = _prepare_metadata(alignment, base_metadata, payload, fasta_bytes)
    metadata_bytes = _encode_metadata(metadata, add_trailing_newline=False)
    version_tuple = _metadata_version_tuple(metadata)
    return EncodedArchive(payload=payload, metadata_bytes=metadata_bytes, metadata=metadata, version_tuple=version_tuple)


def _alignment_to_fasta_bytes(alignment: MultipleSeqAlignment) -> bytes:
    lines: list[str] = []
    for record in alignment:
        identifier = record.id or record.name or "sequence"
        lines.append(f">{identifier}")
        lines.append(str(record.seq))
    return "\n".join(lines).encode("utf-8")


def _prepare_metadata(
    alignment: MultipleSeqAlignment,
    base_metadata: dict[str, object] | None,
    payload: bytes,
    fasta_bytes: bytes,
) -> dict[str, object]:
    metadata = deepcopy(base_metadata) if base_metadata else {}
    metadata.pop("sequence_permutation", None)
    metadata.pop("fallback", None)

    sanitized_keys = [
        "run_length_blocks",
        "max_run_length",
        "columns_with_deviations",
        "bitmask_bytes",
        "bits_per_symbol",
        "payload_encoding",
        "payload_encoded_bytes",
        "payload_raw_bytes",
    ]
    for key in sanitized_keys:
        metadata.pop(key, None)

    ids = [record.id or record.name or f"seq_{index}" for index, record in enumerate(alignment)]
    sequences = [str(record.seq) for record in alignment]
    alphabet = sorted({char for seq in sequences for char in seq})

    format_version = metadata.get("format_version", "0.2.0")
    if isinstance(format_version, (tuple, list)):
        format_version = ".".join(str(part) for part in format_version)
    metadata["format_version"] = format_version
    metadata["codec"] = "ecomp"
    metadata["num_sequences"] = len(sequences)
    metadata["alignment_length"] = alignment.get_alignment_length()
    metadata["alphabet"] = alphabet
    metadata["source_format"] = alignment.annotations.get(
        "source_format", metadata.get("source_format", "fasta")
    )
    metadata["checksum_sha256"] = _checksum(sequences)
    metadata["sequence_ids"] = ids
    metadata["sequence_id_codec"] = "inline"
    metadata["ordering_strategy"] = metadata.get("ordering_strategy", "clipkit")
    metadata["fallback"] = {
        "type": "gzip",
        "format": metadata.get("source_format", "fasta"),
    }
    metadata["payload_encoding"] = "gzip"
    metadata["payload_encoded_bytes"] = len(payload)
    metadata["payload_raw_bytes"] = len(fasta_bytes)

    return metadata


def _metadata_version_tuple(metadata: dict[str, object]) -> tuple[int, int, int]:
    value = metadata.get("format_version", "0.2.0")
    if isinstance(value, str):
        parts = value.split(".")
    elif isinstance(value, Iterable):
        parts = list(value)
    else:
        parts = []

    try:
        major = int(parts[0])
        minor = int(parts[1]) if len(parts) > 1 else 0
        patch = int(parts[2]) if len(parts) > 2 else 0
    except (ValueError, IndexError):
        major, minor, patch = INLINE_METADATA_VERSION
    return major, minor, patch


def _encode_metadata(metadata: dict[str, object], *, add_trailing_newline: bool) -> bytes:
    json_bytes = json.dumps(metadata, sort_keys=True, separators=(",", ":")).encode("utf-8")
    compressed = zlib.compress(json_bytes, level=9)
    if len(compressed) + len(METADATA_COMPRESSED_MAGIC) + 1 < len(json_bytes):
        return METADATA_COMPRESSED_MAGIC + bytes([METADATA_CODEC_VERSION]) + compressed
    if add_trailing_newline:
        return json_bytes + b"\n"
    return json_bytes


def _checksum(sequences: Iterable[str]) -> str:
    digest = sha256()
    for seq in sequences:
        digest.update(seq.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import math
from typing import Union
from Bio.Phylo.BaseTree import Tree

from .ambiguity import (
    ambiguity_map,
    ambiguity_symbols,
    is_nucleotide_alphabet,
    normalize_sequence_type,
    resolved_symbols,
)
from .modes import AmbiguityHandling, StopCodonMode, TrimmingMode

# Always import the standard version for compatibility
from .site_classification import (
    SiteClassificationType,
    determine_site_classification_type,
)

from .settings import DEFAULT_AA_GAP_CHARS
from .stats import TrimmingStats
from .stop_codons import STOP_CODONS, StopCodonMaskingStats

# Counting a batch at a time bounds the temporary integer matrix while still
# letting np.bincount process many columns in one compiled operation.
_COLUMN_COUNT_BATCH_SIZE = 256
_HETEROTACHY_SITE_BATCH_SIZE = 1024
_MAX_DENSE_CODE_SPAN = 512
_MIN_CELLS_FOR_COLUMN_COUNTING = 1_000_000


def _uppercase_records(records: np.ndarray) -> np.ndarray:
    if records.size < 4_096 or records.dtype != np.dtype("U1"):
        return np.char.upper(records)

    code_points = records.view(np.uint32)
    if code_points.max() > 127:
        return np.char.upper(records)

    # Copy the matrix once and bound the temporary masks for large alignments.
    upper = code_points.copy().reshape(-1)
    for start in range(0, upper.size, 1_000_000):
        batch = upper[start : start + 1_000_000]
        lowercase = (batch >= ord("a")) & (batch <= ord("z"))
        np.subtract(batch, 32, out=batch, where=lowercase)
    return upper.view("U1").reshape(records.shape)


def _mostly_dashes(code_points: np.ndarray) -> bool:
    rows, columns = code_points.shape
    sample = code_points[::max(1, rows // 32), ::max(1, columns // 64)]
    return np.count_nonzero(sample != ord("-")) <= sample.size // 8


def _column_character_counts(seq_array: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted character states and their counts in every column."""
    n_rows, n_cols = seq_array.shape
    if n_rows == 0 or n_cols == 0:
        return np.empty(0, dtype="U1"), np.empty((0, n_cols), dtype=np.int64)

    if n_rows <= np.iinfo(np.uint8).max:
        count_dtype = np.uint8
    elif n_rows <= np.iinfo(np.uint16).max:
        count_dtype = np.uint16
    elif n_rows <= np.iinfo(np.uint32).max:
        count_dtype = np.uint32
    else:
        count_dtype = np.uint64

    code_points = seq_array.view(np.uint32).reshape(seq_array.shape)
    min_code = int(code_points.min())
    max_code = int(code_points.max())
    code_span = max_code - min_code + 1

    if code_span <= _MAX_DENSE_CODE_SPAN:
        counts = np.empty((code_span, n_cols), dtype=count_dtype)
        offsets = (
            np.arange(_COLUMN_COUNT_BATCH_SIZE, dtype=np.int64) * code_span - min_code
        )

        gap_code = ord("-")
        sparse_gaps = (
            code_points.size >= 1_000_000
            and min_code <= gap_code <= max_code
            and _mostly_dashes(code_points)
        )

        for start_idx in range(0, n_cols, _COLUMN_COUNT_BATCH_SIZE):
            end_idx = min(start_idx + _COLUMN_COUNT_BATCH_SIZE, n_cols)
            segment = code_points[:, start_idx:end_idx]
            width = segment.shape[1]
            non_gaps = segment != gap_code if sparse_gaps else None
            sparse_batch = non_gaps is not None and (
                np.count_nonzero(non_gaps) <= non_gaps.size // 8
            )
            if sparse_batch:
                rows, columns = np.divmod(np.flatnonzero(non_gaps), width)
                encoded = segment[rows, columns] + offsets[columns]
            else:
                encoded = np.add(segment, offsets[:width], dtype=np.int64)

            histogram = np.bincount(encoded.ravel(), minlength=width * code_span)
            histogram = histogram.reshape(width, code_span)
            if sparse_batch:
                # Account for the omitted dash codes without assigning them to bins.
                histogram[:, gap_code - min_code] = n_rows - histogram.sum(axis=1)
            counts[:, start_idx:end_idx] = histogram.T

        present = np.any(counts, axis=1)
        states = np.array(
            [chr(min_code + offset) for offset in np.flatnonzero(present)],
            dtype="U1",
        )
        return states, counts[present]

    # Unusually wide Unicode alphabets would make the dense code-point table
    # wasteful. Compact them to their observed, sorted states instead.
    states = np.unique(seq_array)
    counts = np.empty((len(states), n_cols), dtype=count_dtype)
    offsets = np.arange(_COLUMN_COUNT_BATCH_SIZE, dtype=np.int64) * len(states)
    for start_idx in range(0, n_cols, _COLUMN_COUNT_BATCH_SIZE):
        end_idx = min(start_idx + _COLUMN_COUNT_BATCH_SIZE, n_cols)
        width = end_idx - start_idx
        encoded = np.searchsorted(states, seq_array[:, start_idx:end_idx])
        encoded += offsets[:width]
        counts[:, start_idx:end_idx] = (
            np.bincount(encoded.ravel(), minlength=width * len(states))
            .reshape(width, len(states))
            .T
        )
    return states, counts


class MSA:
    def __init__(
        self,
        header_info,
        seq_records,
        gap_chars=DEFAULT_AA_GAP_CHARS,
        threads=1,
        requires_uppercase_normalization=False,
        sequence_type="aa",
        ambiguity_handling=AmbiguityHandling.missing,
    ) -> None:
        self.header_info = header_info
        self.seq_records = seq_records
        self._original_length = len(self.seq_records[0])
        self._site_positions_to_keep = np.arange(self._original_length)
        self._site_positions_to_trim = np.array([])
        self._site_classification_types = None
        self._column_character_frequencies = None
        self._column_character_count_cache = {}
        self._observed_sequence_states_cache = None
        self._gap_chars = gap_chars or DEFAULT_AA_GAP_CHARS
        self._codon_size = 3
        self._threads = threads
        self._requires_uppercase_normalization = requires_uppercase_normalization
        self._sequence_type = normalize_sequence_type(sequence_type)
        self._ambiguity_handling = AmbiguityHandling(ambiguity_handling)
        # Cache for expensive computations
        self._site_gappyness_cache = None
        self._site_gap_fraction_cache = None
        self._site_ambiguity_cache = None
        self._site_resolved_fraction_cache = None
        self._site_analyzable_fraction_cache = None
        self._site_entropy_cache = None
        self._site_composition_bias_cache = None

    @staticmethod
    def from_bio_msa(
        alignment: MultipleSeqAlignment,
        gap_chars=None,
        threads=1,
        sequence_type=None,
        ambiguity_handling=AmbiguityHandling.missing,
    ) -> "MSA":
        packed_ascii = (
            type(alignment) is MultipleSeqAlignment
            and len(alignment) > 0
            and len(alignment[0]) >= 4096
            and len(alignment) * len(alignment[0]) >= 1_000_000
            and all(type(record.seq) is Seq for record in alignment)
        )
        header_info = []
        sequences = []
        requires_uppercase_normalization = False
        for rec in alignment:
            header_info.append(
                {"id": rec.id, "name": rec.name, "description": rec.description}
            )
            seq = bytes(rec.seq) if packed_ascii else str(rec.seq)
            if packed_ascii and not seq.isascii():
                # Keep Bio.Seq's ASCII validation for malformed records.
                seq.decode("ascii")
            sequences.append(seq)
            if not requires_uppercase_normalization and seq != seq.upper():
                requires_uppercase_normalization = True

        if sequence_type is None:
            symbols = (
                (map(chr, set(seq)) for seq in sequences)
                if packed_ascii
                else (set(seq) for seq in sequences)
            )
            sequence_type = (
                "nt" if all(is_nucleotide_alphabet(states) for states in symbols) else "aa"
            )

        alignment_length = len(sequences[0])
        if alignment_length == 0:
            seq_records = np.empty((len(sequences), 0), dtype="U1")
        elif packed_ascii:
            code_points = np.empty((len(sequences), alignment_length), dtype=np.uint32)
            for row, sequence in zip(code_points, sequences):
                width = min(len(sequence), alignment_length)
                row[:width] = np.frombuffer(sequence, dtype=np.uint8, count=width)
                if width < alignment_length:
                    row[width:] = 0
            seq_records = code_points.view("U1")
        else:
            seq_records = (
                np.array(sequences, dtype=f"U{alignment_length}")
                .view("U1")
                .reshape(len(sequences), alignment_length)
            )
        return MSA(
            header_info,
            seq_records,
            gap_chars,
            threads,
            requires_uppercase_normalization=requires_uppercase_normalization,
            sequence_type=sequence_type,
            ambiguity_handling=ambiguity_handling,
        )

    def to_bio_msa(self) -> MultipleSeqAlignment:
        sites = (
            self.seq_records
            if isinstance(self.seq_records, np.ndarray)
            and self.length == self._original_length == self.seq_records.shape[1]
            else self.sites_kept
        )
        return self._to_bio_msa(sites)

    def complement_to_bio_msa(self) -> MultipleSeqAlignment:
        return self._to_bio_msa(self.sites_trimmed)

    def mask_stop_codons(self, mode: StopCodonMode) -> StopCodonMaskingStats:
        """Replace selected in-frame stop codons with gap characters."""
        if not isinstance(mode, StopCodonMode):
            mode = StopCodonMode(mode)

        if self._original_length % self._codon_size != 0:
            raise ValueError(
                "Stop codon masking requires an alignment length divisible by 3."
            )
        if self._original_length == 0:
            return StopCodonMaskingStats(mode=mode)

        if "-" not in self._gap_chars:
            self._gap_chars = [*self._gap_chars, "-"]

        normalized_records = (
            _uppercase_records(self.seq_records)
            if self._requires_uppercase_normalization
            else self.seq_records
        )
        codons = normalized_records.reshape(
            self.seq_records.shape[0], -1, self._codon_size
        )
        unavailable_chars = self._effective_unavailable_chars
        gap_chars = np.asarray(sorted(unavailable_chars), dtype="U1")
        if codons.dtype == np.dtype("U1"):
            code_points = codons.view(np.uint32)
            gap_codes = gap_chars.view(np.uint32)
            complete_codons = np.empty(codons.shape[:2], dtype=bool)
            # Bound the temporary arrays used by integer membership checks.
            rows_per_batch = max(1, 262_144 // self._original_length)
            for start in range(0, len(codons), rows_per_batch):
                stop = start + rows_per_batch
                complete_codons[start:stop] = ~np.any(
                    np.isin(code_points[start:stop], gap_codes), axis=2
                )
        else:
            complete_codons = ~np.any(np.isin(codons, gap_chars), axis=2)
        codon_strings = (
            np.ascontiguousarray(codons).view("U3").reshape(codons.shape[:2])
        )
        stop_codons = np.isin(codon_strings, tuple(STOP_CODONS)) & complete_codons

        rows = np.arange(len(codons))
        terminal_positions = complete_codons.shape[1] - 1 - np.argmax(
            complete_codons[:, ::-1], axis=1
        )
        # Rows without complete codons have no stop codons at any position.
        terminal_stops = stop_codons[rows, terminal_positions]
        terminal_masked = 0
        internal_masked = 0
        if mode == StopCodonMode.terminal:
            selected = np.zeros_like(stop_codons)
            selected[rows, terminal_positions] = terminal_stops
            terminal_masked = np.count_nonzero(terminal_stops)
        else:
            selected = stop_codons
            if mode == StopCodonMode.internal:
                selected[rows, terminal_positions] = False
            else:
                terminal_masked = np.count_nonzero(terminal_stops)
            internal_masked = np.count_nonzero(selected) - terminal_masked

        selected_sites = np.repeat(selected, self._codon_size, axis=1)
        self.seq_records[selected_sites] = "-"
        self._reset_analysis_caches()

        return StopCodonMaskingStats(
            mode=mode,
            terminal_masked=int(terminal_masked),
            internal_masked=int(internal_masked),
        )

    def _reset_analysis_caches(self) -> None:
        self._site_classification_types = None
        self._column_character_frequencies = None
        self._column_character_count_cache = {}
        self._observed_sequence_states_cache = None
        self._site_gappyness_cache = None
        self._site_gap_fraction_cache = None
        self._site_ambiguity_cache = None
        self._site_resolved_fraction_cache = None
        self._site_analyzable_fraction_cache = None
        self._site_entropy_cache = None
        self._site_composition_bias_cache = None

    def _to_bio_msa(self, sites) -> MultipleSeqAlignment:
        if sites.shape[1] == 0:
            sequence_rows = [""] * sites.shape[0]
        elif (
            sites.dtype.kind == "U" and sites.dtype.itemsize == np.dtype("U1").itemsize
        ):
            contiguous_sites = np.ascontiguousarray(sites)
            sequence_rows = (
                contiguous_sites.view(f"U{sites.shape[1]}").reshape(-1).tolist()
            )
        else:
            sequence_rows = ["".join(rec) for rec in sites.tolist()]

        # Preserve the full description in the output identifier.
        return MultipleSeqAlignment(
            SeqRecord(Seq(sequence), id=str(info["description"]), description="")
            for sequence, info in zip(sequence_rows, self.header_info)
        )

    def _selected_fasta_rows(self, positions):
        matrix = self.seq_records
        if (
            type(positions) is not np.ndarray
            or positions.ndim != 1
            or len(positions) < 2048
            or not np.can_cast(positions.dtype, np.intp, casting="safe")
            or type(matrix) is not np.ndarray
            or matrix.ndim != 2
            or matrix.dtype != np.dtype("U1")
            or not matrix.size
        ):
            return None
        width = matrix.shape[1]
        code_points = matrix.view(np.uint32)
        # Check before opening the destination so invalid records use the
        # existing materialization path and fail before altering the file.
        if (
            int(positions.min()) < -width
            or int(positions.max()) >= width
            or not all(
                type(info) is dict and type(info.get("description")) is str
                for info in self.header_info
            )
            or code_points.max() > 127
            or code_points.min() == 0
        ):
            return None

        rows = (
            code_points
            if positions is self._site_positions_to_keep
            and len(positions) == self._original_length == width
            else (np.take(row, positions) for row in code_points)
        )
        return (
            (info["description"], row.astype(np.uint8).tobytes())
            for row, info in zip(rows, self.header_info)
        )

    @property
    def trimmed(self):
        if len(self._site_positions_to_trim) == 0:
            return self.seq_records
        return np.delete(self.seq_records, self._site_positions_to_trim, axis=1)

    @property
    def sites_kept(self):
        return np.take(self.seq_records, self._site_positions_to_keep, axis=1)

    @property
    def sites_trimmed(self):
        return np.take(self.seq_records, self._site_positions_to_trim, axis=1)

    @property
    def length(self) -> int:
        return len(self._site_positions_to_keep)

    @property
    def original_length(self):
        return self._original_length

    @property
    def gap_chars(self):
        return self._gap_chars

    @property
    def sequence_type(self) -> str:
        return self._sequence_type

    @property
    def ambiguity_handling(self) -> AmbiguityHandling:
        return self._ambiguity_handling

    @property
    def _normalized_gap_chars(self) -> frozenset[str]:
        return frozenset(str(char).upper() for char in self._gap_chars)

    @property
    def _ambiguity_chars(self) -> frozenset[str]:
        return ambiguity_symbols(self._sequence_type)

    @property
    def _effective_unavailable_chars(self) -> frozenset[str]:
        if self._ambiguity_handling == AmbiguityHandling.missing:
            return self._normalized_gap_chars | self._ambiguity_chars
        return self._normalized_gap_chars

    @property
    def _observed_sequence_states(self) -> np.ndarray:
        if self._observed_sequence_states_cache is None:
            seq_array = (
                _uppercase_records(self.seq_records)
                if self._requires_uppercase_normalization
                else self.seq_records
            )
            self._observed_sequence_states_cache = np.unique(seq_array)
        return self._observed_sequence_states_cache

    def _site_fraction_for_states(
        self, states_to_count: frozenset[str], cache_key: str
    ) -> np.ndarray:
        shared_counts = self._column_character_count_cache.get(
            ("frequencies", self._requires_uppercase_normalization)
        )
        small = self.seq_records.size < _MIN_CELLS_FOR_COLUMN_COUNTING
        direct = (
            shared_counts is None
            and self.seq_records.size >= 1024
            and (small or not self._requires_uppercase_normalization)
            and self.seq_records.dtype == np.dtype("U1")
            and 0 < len(states_to_count) <= (16 if small else 8)
            and all(
                len(state) == 1 and state.isascii() and state != "\0"
                for state in states_to_count
            )
        )
        if shared_counts is None and (small or direct):
            seq_array = (
                _uppercase_records(self.seq_records)
                if self._requires_uppercase_normalization
                else self.seq_records
            )
            if not direct:
                return np.around(
                    np.isin(seq_array, tuple(states_to_count)).mean(axis=0),
                    decimals=4,
                )
            # Count membership directly, with bounded scratch space per row batch.
            code_points = seq_array.view(np.uint32)
            matching_counts = np.zeros(code_points.shape[1], dtype=np.int64)
            rows_per_batch = max(1, 1_000_000 // code_points.shape[1])
            dash = ord("-")
            codes = [ord(state) for state in states_to_count]
            sparse = (
                not small
                and code_points.flags.c_contiguous
                and dash in codes
                and _mostly_dashes(code_points)
            )
            if sparse:
                other_gaps = [code for code in codes if code != dash]
            else:
                low, high = int(code_points.min()), int(code_points.max())
                codes = [code for code in codes if low <= code <= high]
            for start in range(0, len(code_points), rows_per_batch):
                batch = code_points[start : start + rows_per_batch]
                if sparse:
                    non_dash = batch.ravel() != dash
                    if np.count_nonzero(non_dash) <= non_dash.size // 16:
                        indices = np.flatnonzero(non_dash)
                        resolved = indices[~np.isin(batch.ravel()[indices], other_gaps)]
                        # Every cell outside the resolved indices is unavailable.
                        matching_counts += len(batch) - np.bincount(
                            resolved % batch.shape[1], minlength=batch.shape[1]
                        )
                        continue
                matches = np.zeros(batch.shape, dtype=bool)
                for code in codes:
                    matches |= batch == code
                matching_counts += np.count_nonzero(matches, axis=0)
            return np.around(matching_counts / len(code_points), decimals=4)
        if shared_counts is None:
            states, counts = self._get_column_character_counts(
                normalize_case=True, cache_key=cache_key
            )
        else:
            states, counts = shared_counts
        matching_counts = counts[np.isin(states, tuple(states_to_count))].sum(axis=0)
        fractions = matching_counts / self.seq_records.shape[0]
        self._column_character_count_cache.pop(
            (cache_key, self._requires_uppercase_normalization), None
        )
        return np.around(fractions, decimals=4)

    @property
    def site_gap_fraction(self) -> np.ndarray:
        """Fraction of configured gap characters at each site."""
        if self._site_gap_fraction_cache is None:
            self._site_gap_fraction_cache = self._site_fraction_for_states(
                self._normalized_gap_chars, "gap_fraction"
            )
        return self._site_gap_fraction_cache

    @property
    def site_ambiguity(self) -> np.ndarray:
        """Fraction of recognized IUPAC ambiguity symbols at each site."""
        if self._site_ambiguity_cache is None:
            self._site_ambiguity_cache = self._site_fraction_for_states(
                self._ambiguity_chars, "ambiguity_fraction"
            )
        return self._site_ambiguity_cache

    @property
    def site_resolved_fraction(self) -> np.ndarray:
        """Fraction of unambiguous biological states at each site."""
        if self._site_resolved_fraction_cache is None:
            self._site_resolved_fraction_cache = self._site_fraction_for_states(
                resolved_symbols(self._sequence_type), "resolved_fraction"
            )
        return self._site_resolved_fraction_cache

    @property
    def site_analyzable_fraction(self) -> np.ndarray:
        """Fraction of entries contributing to entropy/composition analysis."""
        if self._site_analyzable_fraction_cache is None:
            _, counts = self._get_non_gap_character_counts(
                cache_key="analyzable_fraction"
            )
            self._site_analyzable_fraction_cache = np.around(
                counts.sum(axis=0) / self.seq_records.shape[0], decimals=4
            )
            self._column_character_count_cache.pop(
                (
                    "analyzable_fraction",
                    self._requires_uppercase_normalization,
                ),
                None,
            )
        return self._site_analyzable_fraction_cache

    @property
    def site_gappyness(self) -> np.floating:
        if self._site_gappyness_cache is not None:
            return self._site_gappyness_cache

        self._site_gappyness_cache = self._site_fraction_for_states(
            self._effective_unavailable_chars, "gappyness"
        )
        return self._site_gappyness_cache

    @property
    def is_empty(self) -> bool:
        return np.all(np.take(self.seq_records[0], self._site_positions_to_keep) == "")

    @property
    def site_entropy(self) -> np.ndarray:
        """
        Normalized Shannon entropy per site, excluding gap characters.
        Values are in [0, 1], where 0 means no diversity.
        """
        if self._site_entropy_cache is not None:
            return self._site_entropy_cache

        _, counts_by_state = self._get_non_gap_character_counts(cache_key="entropy")
        total_counts = counts_by_state.sum(axis=0).astype(float)
        if (
            self._ambiguity_handling == AmbiguityHandling.missing
            and self._site_analyzable_fraction_cache is None
        ):
            self._site_analyzable_fraction_cache = np.around(
                total_counts / self.seq_records.shape[0], decimals=4
            )
        observed_states = np.count_nonzero(counts_by_state, axis=0)
        raw_entropy = np.zeros(counts_by_state.shape[1], dtype=float)

        # Iterate over the small alphabet dimension while NumPy evaluates all
        # alignment columns together. This preserves the original summation
        # order without paying Python overhead once per column.
        for state_counts in counts_by_state:
            present = state_counts > 0
            probabilities = state_counts[present] / total_counts[present]
            raw_entropy[present] -= probabilities * np.log2(probabilities)

        entropies = np.zeros(counts_by_state.shape[1], dtype=float)
        variable_sites = observed_states > 1
        entropies[variable_sites] = raw_entropy[variable_sites] / np.log2(
            observed_states[variable_sites]
        )

        self._site_entropy_cache = np.around(entropies, decimals=4)
        self._column_character_count_cache.pop(
            ("entropy", self._requires_uppercase_normalization), None
        )
        return self._site_entropy_cache

    @property
    def site_composition_bias(self) -> np.ndarray:
        """
        Per-site composition bias score in [0, 1], excluding gap characters.
        Higher values indicate stronger dominance by one/few character states.
        """
        if self._site_composition_bias_cache is not None:
            return self._site_composition_bias_cache

        _, counts_by_state = self._get_non_gap_character_counts()
        total_counts = counts_by_state.sum(axis=0).astype(float)
        observed_states = np.count_nonzero(counts_by_state, axis=0)
        dominance = np.zeros(counts_by_state.shape[1], dtype=float)

        # Group columns by the number of observed states so each compact row
        # has the same values and reduction length as the original per-column
        # calculation. This keeps scores identical at four-decimal boundaries.
        for state_count in np.unique(observed_states[observed_states > 1]):
            columns = np.flatnonzero(observed_states == state_count)
            column_counts = counts_by_state[:, columns].T
            present_counts = column_counts[column_counts > 0].reshape(
                len(columns), state_count
            )
            probabilities = present_counts.astype(float) / total_counts[columns, None]
            dominance[columns] = np.sum(np.square(probabilities), axis=1)

        bias_scores = np.zeros(counts_by_state.shape[1], dtype=float)
        bias_scores[observed_states == 1] = 1.0
        variable_sites = observed_states > 1
        min_dominance = 1.0 / observed_states[variable_sites]
        bias_scores[variable_sites] = (dominance[variable_sites] - min_dominance) / (
            1.0 - min_dominance
        )

        self._site_composition_bias_cache = np.around(bias_scores, decimals=4)
        return self._site_composition_bias_cache

    @property
    def stats(self) -> TrimmingStats:
        return TrimmingStats(self)

    def is_any_entry_sequence_only_gaps(self) -> tuple[bool, Union[str, None]]:
        unavailable = tuple(self._effective_unavailable_chars)
        trimmed = self.trimmed
        ascii_gaps = None
        if (
            type(trimmed) is np.ndarray
            and trimmed.ndim == 2
            and trimmed.dtype == np.dtype("U1")
            and trimmed.shape[1] >= 128
            and all(len(char) == 1 and 0 < ord(char) < 128 for char in unavailable)
        ):
            gap_symbols = "".join(unavailable)
            ascii_gaps = (gap_symbols + gap_symbols.lower()).encode("ascii")

        for idx, row in enumerate(trimmed):
            if ascii_gaps is not None and row.view(np.uint32).max() <= 127:
                residues = row.view(np.uint32).astype(np.uint8).tobytes()
                gap_only = not residues.translate(None, ascii_gaps)
            else:
                gap_only = np.all(np.isin(_uppercase_records(row), unavailable))
            if gap_only:
                return True, self.header_info[idx].get("id")
        return False, None

    def trim(
        self,
        mode: TrimmingMode = TrimmingMode.smart_gap,
        gap_threshold=None,
        site_positions_to_trim=None,
        guide_tree: Union[Tree, None] = None,
        codon=False,
        ends_only=False,
    ) -> np.array:
        if site_positions_to_trim is not None:
            if isinstance(site_positions_to_trim, list):
                site_positions_to_trim = np.array(site_positions_to_trim)
            if not isinstance(site_positions_to_trim, np.ndarray):
                raise ValueError("site_positions_to_trim must be a list or np array")

            self._site_positions_to_trim = (
                self.determine_all_codon_sites_to_trim(site_positions_to_trim)
                if codon is True
                else site_positions_to_trim
            )
        else:
            self._site_positions_to_trim = self.determine_site_positions_to_trim(
                mode,
                gap_threshold,
                guide_tree=guide_tree,
                codon=codon,
                ends_only=ends_only,
            )
        if len(self._site_positions_to_trim) > 0:
            self._site_positions_to_keep = np.delete(
                np.arange(self._original_length), self._site_positions_to_trim
            )

    @property
    def column_character_frequencies(self):
        if self._column_character_frequencies is not None:
            return self._column_character_frequencies

        states, counts_by_state = self._get_non_gap_character_counts()
        self._column_character_frequencies = []
        for col_idx in range(counts_by_state.shape[1]):
            present = counts_by_state[:, col_idx] > 0
            self._column_character_frequencies.append(
                dict(
                    zip(
                        states[present],
                        counts_by_state[present, col_idx].tolist(),
                    )
                )
            )

        return self._column_character_frequencies

    def _get_column_character_counts(
        self, normalize_case: bool, cache_key: str = "frequencies"
    ) -> tuple[np.ndarray, np.ndarray]:
        normalize_case = normalize_case and self._requires_uppercase_normalization
        shared_counts = self._column_character_count_cache.get(
            ("frequencies", normalize_case)
        )
        if shared_counts is not None:
            return shared_counts
        key = (cache_key, normalize_case)
        if key not in self._column_character_count_cache:
            seq_array = (
                _uppercase_records(self.seq_records) if normalize_case else self.seq_records
            )
            self._column_character_count_cache[key] = _column_character_counts(
                seq_array
            )
        return self._column_character_count_cache[key]

    def _get_non_gap_character_counts(
        self,
        cache_key: str = "frequencies",
        classification: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        states, counts = self._get_column_character_counts(
            normalize_case=True, cache_key=cache_key
        )
        return self._filter_character_counts(states, counts, classification=classification)

    def _filter_character_counts(
        self, states, counts, *, classification=False, observed_states=None
    ):
        policy = self._ambiguity_handling
        excluded = self._normalized_gap_chars
        if policy == AmbiguityHandling.missing or (
            classification and policy != AmbiguityHandling.literal
        ):
            excluded |= self._ambiguity_chars
        retained = ~np.isin(states, tuple(excluded))
        states = states[retained]
        counts = counts[retained]

        if policy == AmbiguityHandling.fractional and not classification:
            return self._fractionally_expand_ambiguities(
                states, counts, observed_states=observed_states
            )
        return states, counts

    def _fractionally_expand_ambiguities(
        self,
        states: np.ndarray,
        counts: np.ndarray,
        observed_states: np.ndarray = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        expansions = ambiguity_map(
            self._sequence_type,
            states if observed_states is None else observed_states,
        )
        expanded_counts: dict[str, np.ndarray] = {}

        for state, state_counts in zip(states, counts):
            destinations = expansions.get(str(state))
            if destinations is None:
                expanded_counts[str(state)] = (
                    expanded_counts.get(str(state), 0) + state_counts
                )
                continue

            weight = state_counts.astype(float) / len(destinations)
            for destination in destinations:
                expanded_counts[destination] = (
                    expanded_counts.get(destination, 0) + weight
                )

        if not expanded_counts:
            return np.empty(0, dtype="U1"), np.empty((0, counts.shape[1]), dtype=float)

        expanded_states = np.array(sorted(expanded_counts), dtype="U1")
        return expanded_states, np.vstack(
            [expanded_counts[state] for state in expanded_states]
        )

    @property
    def site_classification_types(self):
        if self._site_classification_types is not None:
            return self._site_classification_types

        if self._original_length == 0:
            self._site_classification_types = np.array([])
            return self._site_classification_types

        _, counts_by_state = self._get_non_gap_character_counts(classification=True)
        observed_states = np.count_nonzero(counts_by_state, axis=0)
        states_occurring_at_least_twice = np.count_nonzero(counts_by_state >= 2, axis=0)

        site_classification_types = np.empty(self._original_length, dtype=object)
        site_classification_types.fill(SiteClassificationType.other)
        site_classification_types[
            (states_occurring_at_least_twice == 1) & (observed_states > 1)
        ] = SiteClassificationType.singleton
        site_classification_types[
            (states_occurring_at_least_twice == 1) & (observed_states == 1)
        ] = SiteClassificationType.constant
        site_classification_types[states_occurring_at_least_twice >= 2] = (
            SiteClassificationType.parsimony_informative
        )

        self._site_classification_types = site_classification_types
        return self._site_classification_types

    def determine_site_positions_to_trim(
        self,
        mode,
        gap_threshold,
        guide_tree=None,
        codon=False,
        ends_only=False,
    ):
        if mode in (TrimmingMode.gappy, TrimmingMode.smart_gap):
            sites_to_trim = np.where(self.site_gappyness >= gap_threshold)[0]
        elif mode == TrimmingMode.block_gappy:
            high_gap_sites = np.where(self.site_gappyness >= gap_threshold)[0]
            sites_to_trim = self._retain_contiguous_sites(
                high_gap_sites, min_block_size=2
            )
        elif mode == TrimmingMode.gappyout:
            # Keep sites at the inferred boundary and trim strictly above it.
            sites_to_trim = np.where(self.site_gappyness > gap_threshold)[0]
        elif mode == TrimmingMode.entropy:
            trim_mask = self.site_entropy >= gap_threshold
            if self._ambiguity_handling == AmbiguityHandling.missing:
                trim_mask |= self.site_analyzable_fraction == 0
            sites_to_trim = np.where(trim_mask)[0]
        elif mode == TrimmingMode.composition_bias:
            trim_mask = self.site_composition_bias >= gap_threshold
            if self._ambiguity_handling == AmbiguityHandling.missing:
                trim_mask |= self.site_analyzable_fraction == 0
            sites_to_trim = np.where(trim_mask)[0]
        elif mode == TrimmingMode.heterotachy:
            if guide_tree is None:
                raise ValueError("heterotachy mode requires a guide tree")
            site_heterotachy = self.determine_site_heterotachy(guide_tree)
            trim_mask = site_heterotachy >= gap_threshold
            if self._ambiguity_handling == AmbiguityHandling.missing:
                trim_mask |= self.site_analyzable_fraction == 0
            sites_to_trim = np.where(trim_mask)[0]
        elif mode == TrimmingMode.kpi:
            site_classification_types = self.site_classification_types
            sites_to_trim = np.where(
                site_classification_types
                != SiteClassificationType.parsimony_informative
            )[0]
        elif mode in (TrimmingMode.kpi_gappy, TrimmingMode.kpi_smart_gap):
            site_classification_types = self.site_classification_types
            sites_to_trim = np.flatnonzero(
                (self.site_gappyness > gap_threshold)
                | (
                    site_classification_types
                    != SiteClassificationType.parsimony_informative
                )
            )
        elif mode == TrimmingMode.kpic:
            site_classification_types = self.site_classification_types
            sites_to_trim = np.where(
                (site_classification_types == SiteClassificationType.other)
                | (site_classification_types == SiteClassificationType.singleton)
            )[0]
        elif mode in (TrimmingMode.kpic_gappy, TrimmingMode.kpic_smart_gap):
            site_classification_types = self.site_classification_types
            sites_to_trim = np.flatnonzero(
                (self.site_gappyness >= gap_threshold)
                | (site_classification_types == SiteClassificationType.other)
                | (site_classification_types == SiteClassificationType.singleton)
            )
        elif mode == TrimmingMode.c3:
            sites_to_trim = np.arange(3, self._original_length + 1, 3) - 1
        if codon and mode != TrimmingMode.c3:
            """
            NOTE: ignoring c3 mode otherwise we would ALWAYS trim the entire file by definition.

            For each position in sites_to_trim we need the full triplet of codon positions tuple.
            Example:
                [2, 9] -> [1, 2, 3, 7, 8, 9]
            """
            return self.determine_all_codon_sites_to_trim(sites_to_trim)

        if ends_only:
            sites_to_trim = self.get_consecutive_from_zero_and_max(sites_to_trim)

        return sites_to_trim

    def determine_site_heterotachy(self, guide_tree: Tree) -> np.ndarray:
        """
        Per-site heterotachy score in [0, 1] based on variation in clade-wise
        normalized entropy across internal clades of a guide tree.
        """
        clade_indices = self._get_internal_clade_indices(guide_tree)
        n_sites = self.seq_records.shape[1]
        if not clade_indices:
            return np.zeros(n_sites, dtype=float)

        site_scores = np.zeros(n_sites, dtype=float)
        if self.seq_records.dtype == np.dtype("U1") and len(clade_indices) >= 2:
            observed_states = (
                self._observed_sequence_states
                if self._ambiguity_handling == AmbiguityHandling.fractional
                else None
            )
            for start in range(0, n_sites, _HETEROTACHY_SITE_BATCH_SIZE):
                stop = min(start + _HETEROTACHY_SITE_BATCH_SIZE, n_sites)
                normalized = _uppercase_records(self.seq_records[:, start:stop])
                entropies = np.empty((stop - start, len(clade_indices)), dtype=float)
                for col, indices in enumerate(clade_indices):
                    states, counts = _column_character_counts(normalized[indices])
                    _, counts = self._filter_character_counts(
                        states, counts, observed_states=observed_states
                    )
                    entropies[:, col] = self._entropy_by_observed_state_count(counts)
                site_scores[start:stop] = np.minimum(
                    1.0, np.std(entropies, axis=1) * 2.0
                )
            return np.around(site_scores, decimals=4)

        for site_idx in range(n_sites):
            clade_entropies = [
                self._normalized_entropy_for_clade_site(indices, site_idx)
                for indices in clade_indices
            ]
            if len(clade_entropies) < 2:
                site_scores[site_idx] = 0.0
                continue

            # Standard deviation of values in [0, 1] is bounded by 0.5.
            site_scores[site_idx] = min(1.0, float(np.std(clade_entropies) * 2.0))

        return np.around(site_scores, decimals=4)

    @staticmethod
    def _entropy_by_observed_state_count(counts: np.ndarray) -> np.ndarray:
        # Compact equal-state-count columns so reductions keep scalar order.
        observed = np.count_nonzero(counts, axis=0)
        entropies = np.zeros(counts.shape[1], dtype=float)
        for state_count in np.unique(observed[observed > 1]):
            sites = np.flatnonzero(observed == state_count)
            rows = counts[:, sites].T
            present = rows[rows > 0].reshape(len(sites), state_count)
            probabilities = present.astype(float) / present.sum(axis=1)[:, None]
            raw = -np.sum(probabilities * np.log2(probabilities), axis=1)
            entropies[sites] = raw / math.log2(int(state_count))
        return entropies

    def _get_internal_clade_indices(self, guide_tree: Tree) -> list[list[int]]:
        n_seqs = self.seq_records.shape[0]
        if n_seqs < 3:
            return []

        terminal_to_index = {
            str(info["id"]): idx for idx, info in enumerate(self.header_info)
        }
        clade_indices = []
        seen = set()
        for clade in guide_tree.find_clades(order="level"):
            if clade.is_terminal():
                continue

            taxa_names = [
                str(t.name) for t in clade.get_terminals() if t.name is not None
            ]
            indices = sorted(
                terminal_to_index[name]
                for name in taxa_names
                if name in terminal_to_index
            )
            if len(indices) < 2 or len(indices) >= n_seqs:
                continue

            key = tuple(indices)
            if key in seen:
                continue
            seen.add(key)
            clade_indices.append(indices)
        return clade_indices

    def _normalized_entropy_for_clade_site(
        self,
        clade_indices: list[int],
        site_idx: int,
    ) -> float:
        column = np.char.upper(self.seq_records[clade_indices, site_idx])
        states, counts = np.unique(column, return_counts=True)
        counts = counts.reshape(-1, 1)

        observed_states = (
            self._observed_sequence_states
            if self._ambiguity_handling == AmbiguityHandling.fractional
            else None
        )
        _, counts = self._filter_character_counts(
            states, counts, observed_states=observed_states
        )

        present_counts = counts[:, 0]
        present_counts = present_counts[present_counts > 0]
        if len(present_counts) <= 1:
            return 0.0

        probs = present_counts.astype(float) / float(np.sum(present_counts))
        raw_entropy = -float(np.sum(probs * np.log2(probs)))
        max_entropy = math.log2(len(present_counts))
        return raw_entropy / max_entropy if max_entropy > 0 else 0.0

    @staticmethod
    def _retain_contiguous_sites(
        positions: np.ndarray,
        min_block_size: int = 2,
    ) -> np.ndarray:
        if positions.size == 0:
            return positions

        contiguous_blocks = []
        start = 0

        for idx in range(1, len(positions) + 1):
            is_break = idx == len(positions) or positions[idx] != positions[idx - 1] + 1
            if not is_break:
                continue

            block = positions[start:idx]
            if len(block) >= min_block_size:
                contiguous_blocks.append(block)
            start = idx

        if not contiguous_blocks:
            return np.array([], dtype=positions.dtype)

        return np.concatenate(contiguous_blocks)

    def determine_gappyout_gap_threshold(self) -> float:
        """
        Infer a gap threshold from the empirical gappyness distribution.

        This gappyout-inspired heuristic finds the largest jump between
        sorted per-site gappyness values and sets the threshold to the
        midpoint of that jump. If no clear jump exists, it falls back to
        a conservative high-quantile threshold.
        """
        gappyness = np.sort(self.site_gappyness.astype(float))
        if gappyness.size == 0:
            return 1.0

        # Uniform distributions do not expose a natural split.
        if np.allclose(gappyness, gappyness[0]):
            return 1.0

        jumps = np.diff(gappyness)
        max_jump_idx = int(np.argmax(jumps))
        max_jump = float(jumps[max_jump_idx])

        if max_jump >= 0.1:
            lower = float(gappyness[max_jump_idx])
            upper = float(gappyness[max_jump_idx + 1])
            threshold = (lower + upper) / 2.0
        else:
            threshold = float(np.quantile(gappyness, 0.95))
            threshold = max(threshold, 0.9)

        return round(max(0.0, min(1.0, threshold)), 4)

    def generate_debug_log_info(self):
        """
        Returns tuples of site position, keep or trim, site classification type, and gappyness
        """
        keep_or_trim_lookup = {}
        for keep_idx in self._site_positions_to_keep:
            keep_or_trim_lookup[keep_idx] = "keep"
        for trim_idx in self._site_positions_to_trim:
            keep_or_trim_lookup[trim_idx] = "trim"

        for idx, gappyness in enumerate(self.site_gappyness):
            yield (
                idx,
                keep_or_trim_lookup[idx],
                self.site_classification_types[idx],
                gappyness,
            )

    def determine_all_codon_sites_to_trim(self, sites_to_trim):
        """
        For each position in sites_to_trim we need the full triplet of codon positions.

        Sites to trim -> all codon sites to trim
        [2, 8] -> [0, 1, 2, 6, 7, 8]
        """
        # Vectorized codon triplet calculation
        if len(sites_to_trim) == 0:
            return np.array([])

        sites_to_trim = np.asarray(sites_to_trim)
        blocks = np.unique(sites_to_trim // self._codon_size)
        triplets = (
            blocks[:, None] * self._codon_size + np.arange(self._codon_size)
        ).reshape(-1)
        return triplets[triplets < self._original_length]

    def determine_codon_triplet_positions(self, alignment_position):
        """
        Block 0 -> [0,1,2], block 1 -> [3,4,5]

        We filter to make sure we are not including any positions out of range
        """
        block = alignment_position // self._codon_size
        codon_triplet_index_start = block * self._codon_size
        sites = [
            codon_triplet_index_start,
            codon_triplet_index_start + 1,
            codon_triplet_index_start + 2,
        ]
        return list(filter(lambda x: x <= self._original_length - 1, sites))

    def get_consecutive_starting_from_zero(self, arr):
        if arr.size == 0:
            return arr

        zero_locs = np.where(arr == 0)[0]
        if zero_locs.size == 0:
            return np.array([], dtype=arr.dtype)

        start_idx = zero_locs[0]

        run = [arr[start_idx]]
        i = start_idx
        while i + 1 < len(arr) and arr[i + 1] == arr[i] + 1:
            run.append(arr[i + 1])
            i += 1

        return np.array(run, dtype=arr.dtype)

    def get_consecutive_ending_with_max(self, arr):
        if arr.size == 0:
            return arr

        max_locs = np.where(arr == (self._original_length - 1))[0]
        try:
            end_idx = max_locs[-1]
        except IndexError:
            return np.array([], dtype=arr.dtype)
        run = [arr[end_idx]]
        i = end_idx
        while i - 1 >= 0 and arr[i - 1] == arr[i] - 1:
            run.append(arr[i - 1])
            i -= 1

        run.reverse()

        return np.array(run, dtype=arr.dtype)

    def get_consecutive_from_zero_and_max(self, arr):
        run_from_zero = self.get_consecutive_starting_from_zero(arr)
        run_with_max = self.get_consecutive_ending_with_max(arr)

        return np.concatenate([run_from_zero, run_with_max])

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import math
from typing import Union
from Bio.Phylo.BaseTree import Tree

from .modes import StopCodonMode, TrimmingMode

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
_COLUMN_COUNT_BATCH_SIZE = 1024
_MAX_DENSE_CODE_SPAN = 512
_MIN_CELLS_FOR_COLUMN_COUNTING = 1_000_000


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
        offsets = np.arange(_COLUMN_COUNT_BATCH_SIZE, dtype=np.int64) * code_span

        for start_idx in range(0, n_cols, _COLUMN_COUNT_BATCH_SIZE):
            end_idx = min(start_idx + _COLUMN_COUNT_BATCH_SIZE, n_cols)
            width = end_idx - start_idx
            encoded = code_points[:, start_idx:end_idx].astype(np.int64)
            encoded -= min_code
            encoded += offsets[:width]
            counts[:, start_idx:end_idx] = (
                np.bincount(encoded.ravel(), minlength=width * code_span)
                .reshape(width, code_span)
                .T
            )

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
    ) -> None:
        self.header_info = header_info
        self.seq_records = seq_records
        self._original_length = len(self.seq_records[0])
        self._site_positions_to_keep = np.arange(self._original_length)
        self._site_positions_to_trim = np.array([])
        self._site_classification_types = None
        self._column_character_frequencies = None
        self._column_character_count_cache = {}
        self._gap_chars = gap_chars or DEFAULT_AA_GAP_CHARS
        self._codon_size = 3
        self._threads = threads
        self._requires_uppercase_normalization = requires_uppercase_normalization
        # Cache for expensive computations
        self._site_gappyness_cache = None
        self._site_entropy_cache = None
        self._site_composition_bias_cache = None

    @staticmethod
    def from_bio_msa(
        alignment: MultipleSeqAlignment, gap_chars=None, threads=1
    ) -> "MSA":
        header_info = []
        sequences = []
        requires_uppercase_normalization = False
        for rec in alignment:
            header_info.append(
                {"id": rec.id, "name": rec.name, "description": rec.description}
            )
            seq = str(rec.seq)
            sequences.append(seq)
            if not requires_uppercase_normalization and seq != seq.upper():
                requires_uppercase_normalization = True

        alignment_length = len(sequences[0])
        if alignment_length == 0:
            seq_records = np.empty((len(sequences), 0), dtype="U1")
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
        )

    def to_bio_msa(self) -> MultipleSeqAlignment:
        return self._to_bio_msa(self.sites_kept)

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

        codons = np.char.upper(self.seq_records).reshape(
            self.seq_records.shape[0], -1, self._codon_size
        )
        gap_chars = np.char.upper(np.asarray(self._gap_chars, dtype="U1"))
        complete_codons = ~np.any(np.isin(codons, gap_chars), axis=2)
        codon_strings = (
            np.ascontiguousarray(codons).view("U3").reshape(codons.shape[:2])
        )
        stop_codons = np.isin(codon_strings, tuple(STOP_CODONS)) & complete_codons

        terminal_to_mask = np.zeros(stop_codons.shape, dtype=bool)
        internal_to_mask = np.zeros(stop_codons.shape, dtype=bool)

        for sequence_index in range(self.seq_records.shape[0]):
            complete_positions = np.flatnonzero(complete_codons[sequence_index])
            if complete_positions.size == 0:
                continue

            terminal_position = complete_positions[-1]
            terminal_to_mask[sequence_index, terminal_position] = stop_codons[
                sequence_index, terminal_position
            ]
            internal_to_mask[sequence_index] = stop_codons[sequence_index]
            internal_to_mask[sequence_index, terminal_position] = False

        selected = np.zeros(stop_codons.shape, dtype=bool)
        if mode in (StopCodonMode.terminal, StopCodonMode.all):
            selected |= terminal_to_mask
        if mode in (StopCodonMode.internal, StopCodonMode.all):
            selected |= internal_to_mask

        selected_sites = np.repeat(selected, self._codon_size, axis=1)
        self.seq_records[selected_sites] = "-"
        self._reset_analysis_caches()

        return StopCodonMaskingStats(
            mode=mode,
            terminal_masked=int(np.count_nonzero(terminal_to_mask & selected)),
            internal_masked=int(np.count_nonzero(internal_to_mask & selected)),
        )

    def _reset_analysis_caches(self) -> None:
        self._site_classification_types = None
        self._column_character_frequencies = None
        self._column_character_count_cache = {}
        self._site_gappyness_cache = None
        self._site_entropy_cache = None
        self._site_composition_bias_cache = None

    def _to_bio_msa(self, sites) -> MultipleSeqAlignment:
        # NOTE: we use the description as the id to preserve the full sequence description - see issue #20
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

        return MultipleSeqAlignment(
            [
                SeqRecord(Seq(rec), id=str(info["description"]), description="")
                for rec, info in zip(sequence_rows, self.header_info)
            ]
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
    def site_gappyness(self) -> np.floating:
        if self._site_gappyness_cache is not None:
            return self._site_gappyness_cache

        if self.seq_records.size < _MIN_CELLS_FOR_COLUMN_COUNTING:
            site_gappyness = np.isin(self.seq_records, self._gap_chars).mean(axis=0)
        else:
            states, counts = self._get_column_character_counts(
                normalize_case=False, cache_key="gappyness"
            )
            gap_counts = counts[np.isin(states, self._gap_chars)].sum(axis=0)
            site_gappyness = gap_counts / self.seq_records.shape[0]
        self._site_gappyness_cache = np.around(site_gappyness, decimals=4)
        self._column_character_count_cache.pop(("gappyness", False), None)
        return self._site_gappyness_cache

    @property
    def is_empty(self) -> bool:
        all_zeros = np.all(self.sites_kept[0] == "")
        return all_zeros

    @property
    def site_entropy(self) -> np.ndarray:
        """
        Normalized Shannon entropy per site, excluding gap characters.
        Values are in [0, 1], where 0 means no diversity.
        """
        if self._site_entropy_cache is not None:
            return self._site_entropy_cache

        _, counts_by_state = self._get_non_gap_character_counts(
            normalize_gap_chars=True, cache_key="entropy"
        )
        n_cols = counts_by_state.shape[1]
        entropies = np.zeros(n_cols, dtype=float)

        for col_idx in range(n_cols):
            non_gap_counts = counts_by_state[:, col_idx]
            non_gap_counts = non_gap_counts[non_gap_counts > 0]

            if len(non_gap_counts) <= 1:
                continue

            total = float(sum(non_gap_counts))
            probs = [count / total for count in non_gap_counts]
            raw_entropy = -sum(p * math.log2(p) for p in probs if p > 0.0)
            max_entropy = math.log2(len(non_gap_counts))
            entropies[col_idx] = raw_entropy / max_entropy if max_entropy > 0 else 0.0

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
        bias_scores = np.zeros(counts_by_state.shape[1], dtype=float)

        for idx in range(counts_by_state.shape[1]):
            counts = counts_by_state[:, idx]
            counts = counts[counts > 0].astype(float)
            total = counts.sum()

            if total == 0:
                bias_scores[idx] = 0.0
                continue

            n_states = len(counts)
            if n_states == 1:
                bias_scores[idx] = 1.0
                continue

            probs = counts / total
            dominance = float(np.sum(np.square(probs)))
            min_dominance = 1.0 / float(n_states)
            bias_scores[idx] = (dominance - min_dominance) / (1.0 - min_dominance)

        self._site_composition_bias_cache = np.around(bias_scores, decimals=4)
        return self._site_composition_bias_cache

    @property
    def stats(self) -> TrimmingStats:
        return TrimmingStats(self)

    def is_any_entry_sequence_only_gaps(self) -> tuple[bool, Union[str, None]]:
        for idx, row in enumerate(self.trimmed):
            if np.all(row == row[0]) and (  # all values the same
                row[0] in self.gap_chars
            ):
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
                        counts_by_state[present, col_idx].astype(np.int64),
                    )
                )
            )

        return self._column_character_frequencies

    def _get_column_character_counts(
        self, normalize_case: bool, cache_key: str = "frequencies"
    ) -> tuple[np.ndarray, np.ndarray]:
        normalize_case = normalize_case and self._requires_uppercase_normalization
        key = (cache_key, normalize_case)
        if key not in self._column_character_count_cache:
            seq_array = (
                np.char.upper(self.seq_records) if normalize_case else self.seq_records
            )
            self._column_character_count_cache[key] = _column_character_counts(
                seq_array
            )
        return self._column_character_count_cache[key]

    def _get_non_gap_character_counts(
        self,
        normalize_gap_chars: bool = False,
        cache_key: str = "frequencies",
    ) -> tuple[np.ndarray, np.ndarray]:
        states, counts = self._get_column_character_counts(
            normalize_case=True, cache_key=cache_key
        )
        gap_chars = (
            [gap.upper() for gap in self._gap_chars]
            if normalize_gap_chars
            else self._gap_chars
        )
        non_gap_states = ~np.isin(states, gap_chars)
        return states[non_gap_states], counts[non_gap_states]

    @property
    def site_classification_types(self):
        if self._site_classification_types is not None:
            return self._site_classification_types

        if self._original_length == 0:
            self._site_classification_types = np.array([])
            return self._site_classification_types

        _, counts_by_state = self._get_non_gap_character_counts()
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
            sites_to_trim = np.where(self.site_entropy >= gap_threshold)[0]
        elif mode == TrimmingMode.composition_bias:
            sites_to_trim = np.where(self.site_composition_bias >= gap_threshold)[0]
        elif mode == TrimmingMode.heterotachy:
            if guide_tree is None:
                raise ValueError("heterotachy mode requires a guide tree")
            site_heterotachy = self.determine_site_heterotachy(guide_tree)
            sites_to_trim = np.where(site_heterotachy >= gap_threshold)[0]
        elif mode == TrimmingMode.kpi:
            site_classification_types = self.site_classification_types
            sites_to_trim = np.where(
                site_classification_types
                != SiteClassificationType.parsimony_informative
            )[0]
        elif mode in (TrimmingMode.kpi_gappy, TrimmingMode.kpi_smart_gap):
            sites_to_trim_gaps_based = np.where(self.site_gappyness > gap_threshold)[0]
            site_classification_types = self.site_classification_types
            sites_to_trim_classification_based = np.where(
                site_classification_types
                != SiteClassificationType.parsimony_informative
            )[0]
            sites_to_trim = np.unique(
                np.concatenate(
                    (sites_to_trim_gaps_based, sites_to_trim_classification_based)
                )
            )
        elif mode == TrimmingMode.kpic:
            site_classification_types = self.site_classification_types
            sites_to_trim = np.where(
                (site_classification_types == SiteClassificationType.other)
                | (site_classification_types == SiteClassificationType.singleton)
            )[0]
        elif mode in (TrimmingMode.kpic_gappy, TrimmingMode.kpic_smart_gap):
            sites_to_trim_gaps_based = np.where(self.site_gappyness >= gap_threshold)[0]

            site_classification_types = self.site_classification_types
            sites_to_trim_classification_based = np.where(
                (site_classification_types == SiteClassificationType.other)
                | (site_classification_types == SiteClassificationType.singleton)
            )[0]

            sites_to_trim = np.unique(
                np.concatenate(
                    (sites_to_trim_gaps_based, sites_to_trim_classification_based)
                )
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
        column = self.seq_records[clade_indices, site_idx]
        gap_chars_upper = {char.upper() for char in self._gap_chars}
        non_gap = [char for char in column if char.upper() not in gap_chars_upper]

        if len(non_gap) <= 1:
            return 0.0

        unique, counts = np.unique(non_gap, return_counts=True)
        if len(unique) <= 1:
            return 0.0

        probs = counts.astype(float) / float(np.sum(counts))
        raw_entropy = -float(np.sum(probs * np.log2(probs)))
        max_entropy = math.log2(len(unique))
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
        blocks = sites_to_trim // self._codon_size

        # Create all triplet positions for each block
        triplets = []
        for block in np.unique(blocks):
            start = block * self._codon_size
            positions = np.arange(
                start, min(start + self._codon_size, self._original_length)
            )
            triplets.extend(positions)

        return np.unique(triplets)

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

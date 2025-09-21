from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from itertools import chain
from typing import Union
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count
from functools import lru_cache

from .modes import TrimmingMode
# Always import the standard version for compatibility
from .site_classification import (
    SiteClassificationType,
    determine_site_classification_type,
)

# Optimizations are built-in, no need for external modules
USE_OPTIMIZED = True
batch_classify_sites_wrapper = None
calculate_column_frequency_fast = None
from .settings import DEFAULT_AA_GAP_CHARS
from .stats import TrimmingStats


# Vectorized helper functions
def _vectorized_column_frequencies(seq_array, gap_chars):
    """
    Vectorized calculation of character frequencies for all columns at once
    using NumPy operations.
    """
    n_cols = seq_array.shape[1]
    column_frequencies = []

    # Convert to uppercase once
    seq_array_upper = np.char.upper(seq_array.astype(str))

    for i in range(n_cols):
        column = seq_array_upper[:, i]
        unique, counts = np.unique(column, return_counts=True)
        freqs = dict(zip(unique, counts))

        # Remove gap characters
        for gap_char in gap_chars:
            freqs.pop(gap_char, None)

        column_frequencies.append(freqs)

    return column_frequencies


def _batch_column_frequencies(seq_array, gap_chars, batch_size=100):
    """
    Process columns in batches for better memory efficiency
    """
    n_cols = seq_array.shape[1]
    column_frequencies = []

    # Convert to uppercase once for the entire array
    seq_array_upper = np.char.upper(seq_array.astype(str))

    for start_idx in range(0, n_cols, batch_size):
        end_idx = min(start_idx + batch_size, n_cols)
        batch = seq_array_upper[:, start_idx:end_idx]

        for i in range(batch.shape[1]):
            column = batch[:, i]
            unique, counts = np.unique(column, return_counts=True)
            freqs = dict(zip(unique, counts))

            # Remove gap characters
            for gap_char in gap_chars:
                freqs.pop(gap_char, None)

            column_frequencies.append(freqs)

    return column_frequencies


# Module-level helper functions for parallel processing
def _calculate_column_frequency_helper(args):
    """Helper function for parallel processing of column frequencies"""
    column, gap_chars = args
    col_sorted_unique_values_for, col_counts_per_char = np.unique(
        [char.upper() for char in column], return_counts=True
    )
    freqs = dict(zip(col_sorted_unique_values_for, col_counts_per_char))
    for gap_char in gap_chars:
        try:
            del freqs[gap_char]
        except KeyError:
            continue
    return freqs


class MSA:
    def __init__(
        self, header_info, seq_records, gap_chars=DEFAULT_AA_GAP_CHARS, threads=1
    ) -> None:
        self.header_info = header_info
        self.seq_records = seq_records
        self._original_length = len(self.seq_records[0])
        self._site_positions_to_keep = np.arange(self._original_length)
        self._site_positions_to_trim = np.array([])
        self._site_classification_types = None
        self._column_character_frequencies = None
        self._gap_chars = gap_chars
        self._codon_size = 3
        self._threads = threads
        # Cache for expensive computations
        self._site_gappyness_cache = None

    @staticmethod
    def from_bio_msa(alignment: MultipleSeqAlignment, gap_chars=None, threads=1) -> "MSA":
        header_info = [
            {"id": rec.id, "name": rec.name, "description": rec.description}
            for rec in alignment
        ]
        seq_records = np.array([list(rec) for rec in alignment])
        return MSA(header_info, seq_records, gap_chars, threads)

    def to_bio_msa(self) -> MultipleSeqAlignment:
        return self._to_bio_msa(self.sites_kept)

    def complement_to_bio_msa(self) -> MultipleSeqAlignment:
        return self._to_bio_msa(self.sites_trimmed)

    def _to_bio_msa(self, sites) -> MultipleSeqAlignment:
        # NOTE: we use the description as the id to preserve the full sequence description - see issue #20
        return MultipleSeqAlignment(
            [
                SeqRecord(
                    Seq("".join(rec)), id=str(info["description"]), description=""
                )
                for rec, info in zip(sites.tolist(), self.header_info)
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

        # Vectorized calculation using broadcasting
        site_gappyness = (np.isin(self.seq_records, self._gap_chars)).mean(axis=0)
        self._site_gappyness_cache = np.around(site_gappyness, decimals=4)
        return self._site_gappyness_cache

    @property
    def is_empty(self) -> bool:
        all_zeros = np.all(self.sites_kept[0] == "")
        return all_zeros

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
                codon,
                ends_only,
            )
        if len(self._site_positions_to_trim) > 0:
            self._site_positions_to_keep = np.delete(
                np.arange(self._original_length), self._site_positions_to_trim
            )

    @property
    def column_character_frequencies(self):
        if self._column_character_frequencies is not None:
            return self._column_character_frequencies

        # Adaptive threshold based on alignment size and thread count
        # Lower threshold for better parallel efficiency
        parallel_threshold = max(1000, 5000 // self._threads)

        # Use vectorized or parallel processing based on size
        if self._original_length < 500:
            # Small alignments: use vectorized numpy operations
            self._column_character_frequencies = _vectorized_column_frequencies(
                self.seq_records, self.gap_chars
            )
        elif self._threads > 1 and self._original_length > parallel_threshold:
            # Large alignments with multiple threads: use parallel processing
            columns = self.seq_records.T

            # Dynamic worker calculation
            n_workers = min(self._threads, max(2, min(self._original_length // 500, cpu_count())))

            # Prepare arguments for parallel processing
            args_list = [(column, self.gap_chars) for column in columns]

            # Use ProcessPoolExecutor with optimized chunk size
            chunksize = max(10, len(args_list) // (n_workers * 10))

            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                self._column_character_frequencies = list(executor.map(
                    _calculate_column_frequency_helper,
                    args_list,
                    chunksize=chunksize
                ))
        else:
            # Medium alignments or single-threaded: use batch processing
            self._column_character_frequencies = _batch_column_frequencies(
                self.seq_records, self.gap_chars, batch_size=200
            )

        # Ensure we always return the frequencies (not None)
        if self._column_character_frequencies is None:
            self._column_character_frequencies = _batch_column_frequencies(
                self.seq_records, self.gap_chars, batch_size=200
            )

        return self._column_character_frequencies

    @property
    def site_classification_types(self):
        if self._site_classification_types is not None:
            return self._site_classification_types

        col_char_freqs = self.column_character_frequencies

        # Adaptive parallel threshold
        parallel_threshold = max(1000, 5000 // self._threads)

        if self._threads > 1 and len(col_char_freqs) > parallel_threshold:
            n_workers = min(self._threads, max(2, min(len(col_char_freqs) // 500, cpu_count())))
            chunksize = max(10, len(col_char_freqs) // (n_workers * 10))

            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                site_classification_types = list(executor.map(
                    determine_site_classification_type,
                    col_char_freqs,
                    chunksize=chunksize
                ))
            site_classification_types = np.array(site_classification_types)
        else:
            # Vectorized single-threaded implementation
            site_classification_types = np.array(
                [determine_site_classification_type(freq) for freq in col_char_freqs]
            )

        self._site_classification_types = site_classification_types
        return self._site_classification_types

    def determine_site_positions_to_trim(
        self,
        mode,
        gap_threshold,
        codon=False,
        ends_only=False,
    ):
        if mode in (TrimmingMode.gappy, TrimmingMode.smart_gap):
            sites_to_trim = np.where(self.site_gappyness >= gap_threshold)[0]
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
            positions = np.arange(start, min(start + self._codon_size, self._original_length))
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
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from itertools import chain
from typing import Union

from .modes import TrimmingMode
from .site_classification import (
    SiteClassificationType,
    determine_site_classification_type,
)
from .settings import DEFAULT_AA_GAP_CHARS
from .stats import TrimmingStats


class MSA:
    def __init__(
        self, header_info, seq_records, gap_chars=DEFAULT_AA_GAP_CHARS
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

    @staticmethod
    def from_bio_msa(alignment: MultipleSeqAlignment, gap_chars=None) -> "MSA":
        header_info = [
            {"id": rec.id, "name": rec.name, "description": rec.description}
            for rec in alignment
        ]
        seq_records = np.array([list(rec) for rec in alignment])
        return MSA(header_info, seq_records, gap_chars)

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
        site_gappyness = (np.isin(self.seq_records, self._gap_chars)).mean(axis=0)
        return np.around(site_gappyness, decimals=4)

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
            )
        self._site_positions_to_keep = np.delete(
            np.arange(self._original_length), self._site_positions_to_trim
        )

    @property
    def column_character_frequencies(self):
        if self._column_character_frequencies is not None:
            return self._column_character_frequencies

        column_character_frequencies = []
        for column in self.seq_records.T:
            col_sorted_unique_values_for, col_counts_per_char = np.unique(
                [char.upper() for char in column], return_counts=True
            )
            freqs = dict(zip(col_sorted_unique_values_for, col_counts_per_char))
            for gap_char in self.gap_chars:
                try:
                    del freqs[gap_char]
                except KeyError:
                    continue
            column_character_frequencies.append(freqs)
        self._column_character_frequencies = column_character_frequencies
        return self._column_character_frequencies

    @property
    def site_classification_types(self):
        if self._site_classification_types is not None:
            return self._site_classification_types

        site_classification_types = np.array(
            [
                determine_site_classification_type(col_char_freq)
                for col_char_freq in self.column_character_frequencies
            ]
        )
        self._site_classification_types = site_classification_types
        return self._site_classification_types

    def determine_site_positions_to_trim(self, mode, gap_threshold, codon=False):
        if mode in (TrimmingMode.gappy, TrimmingMode.smart_gap):
            sites_to_trim = np.where(self.site_gappyness >= gap_threshold)[0]
        elif mode == TrimmingMode.kpi:
            # col_char_freqs = self.column_character_frequencies
            site_classification_types = self.site_classification_types
            sites_to_trim = np.where(
                site_classification_types
                != SiteClassificationType.parsimony_informative
            )[0]
        elif mode in (TrimmingMode.kpi_gappy, TrimmingMode.kpi_smart_gap):
            sites_to_trim_gaps_based = np.where(self.site_gappyness > gap_threshold)[0]
            # col_char_freqs = self.column_character_frequencies
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
            # col_char_freqs = self.column_character_frequencies
            site_classification_types = self.site_classification_types
            sites_to_trim = np.where(
                (site_classification_types == SiteClassificationType.other)
                | (site_classification_types == SiteClassificationType.singleton)
            )[0]
        elif mode in (TrimmingMode.kpic_gappy, TrimmingMode.kpic_smart_gap):
            sites_to_trim_gaps_based = np.where(self.site_gappyness >= gap_threshold)[0]

            # col_char_freqs = self.column_character_frequencies
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

        if codon:
            """
            For each position in sites_to_trim we need the full triplet of codon positions tuple.
            Example:
                [2, 9] -> [1, 2, 3, 7, 8, 9]
            """
            return self.determine_all_codon_sites_to_trim(sites_to_trim)

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
        sites_to_trim_codon = [
            self.determine_codon_triplet_positions(site_pos)
            for site_pos in sites_to_trim
        ]
        flattened_unique_sites = list(set(chain(*sites_to_trim_codon)))
        return np.array(flattened_unique_sites)

    def determine_codon_triplet_positions(self, alignment_position):
        """
        Block 0 -> [0,1,2], block 1 -> [3,4,5]

        We filter to make sure we are not including any positions out of range
        """
        block = alignment_position // 3
        codon_triplet_index_start = block * 3
        sites = [
            codon_triplet_index_start,
            codon_triplet_index_start + 1,
            codon_triplet_index_start + 2,
        ]
        return list(filter(lambda x: x <= self._original_length - 1, sites))

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from typing import Union

from .modes import SiteClassificationType, TrimmingMode
from .site_classification import determine_site_classification_type
from .settings import DEFAULT_AA_GAP_CHARS
from .stats import TrimmingStats


class MSA:
    def __init__(self, header_info, seq_records) -> None:
        self.header_info = header_info
        self.seq_records = seq_records
        self._original_length = len(self.seq_records[0])
        self._site_positions_to_keep = np.arange(self._original_length)
        self._site_positions_to_trim = np.array([])
        self._gap_chars = DEFAULT_AA_GAP_CHARS
        # self.site_classification_counts = None TODO:

    @property
    def trimmed(self):
        return np.delete(self.seq_records, self._site_positions_to_trim, axis=1)[0]

    @property
    def sites_kept(self):
        return np.take(self.seq_records, self._site_positions_to_keep, axis=1)

    @property
    def sites_trimmed(self):
        return np.take(self.seq_records, self._site_positions_to_trim, axis=1)

    @property
    def length(self) -> int:
        return self._original_length if self.trimmed is None else len(self._site_positions_to_keep)

    @property
    def original_length(self):
        return self._original_length

    @property
    def gap_chars(self):
        return self._gap_chars

    @gap_chars.setter
    def gap_chars(self, gap_chars):
        self._gap_chars = gap_chars

    @property
    def site_gappyness(self) -> np.floating:
        return (np.isin(self.seq_records, self._gap_chars)).mean(axis=0)

    @property
    def is_empty(self) -> bool:
        print(self.sites_kept)
        all_zeros = np.all(self.sites_kept[0] == b"")
        return all_zeros

    def trim(
        self,
        mode: TrimmingMode,
        gap_threshold=None,
    ) -> np.array:
        self._site_positions_to_trim = self.determine_site_positions_to_trim(mode, gap_threshold)
        self._site_positions_to_keep = np.delete(np.arange(self._original_length), self._site_positions_to_trim)

    def column_character_frequencies(self):
        column_character_frequencies = []
        for column in self.seq_records.T:
            col_sorted_unique_values_for, col_counts_per_char = np.unique(
                column, return_counts=True
            )
            freqs = dict(zip(col_sorted_unique_values_for, col_counts_per_char))
            for gap_char in self.gap_chars:
                try:
                    del freqs[gap_char]
                except KeyError:
                    continue
            column_character_frequencies.append(freqs)
        return column_character_frequencies

    def determine_site_positions_to_trim(self, mode, gap_threshold):
        if mode in (TrimmingMode.gappy, TrimmingMode.smart_gap):
            sites_to_trim = np.where(self.site_gappyness > gap_threshold)[0]
        elif mode == TrimmingMode.kpi:
            col_char_freqs = self.column_character_frequencies()
            site_classification_types = np.array(
                [
                    determine_site_classification_type(col_char_freq)
                    for col_char_freq in col_char_freqs
                ]
            )
            sites_to_trim = np.where(
                site_classification_types
                != SiteClassificationType.parsimony_informative
            )[0]
        elif mode in (TrimmingMode.kpi_gappy, TrimmingMode.kpi_smart_gap):
            sites_to_trim_gaps_based = np.where(self.site_gappyness > gap_threshold)[0]
            col_char_freqs = self.column_character_frequencies()
            site_classification_types = np.array(
                [
                    determine_site_classification_type(col_char_freq)
                    for col_char_freq in col_char_freqs
                ]
            )
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
            col_char_freqs = self.column_character_frequencies()
            site_classification_types = np.array(
                [
                    determine_site_classification_type(col_char_freq)
                    for col_char_freq in col_char_freqs
                ]
            )
            sites_to_trim = np.where(
                (
                    site_classification_types
                    != SiteClassificationType.parsimony_informative
                )
                | (site_classification_types != SiteClassificationType.constant)
            )[0]
        elif mode in (TrimmingMode.kpic_gappy, TrimmingMode.kpic_smart_gap):
            sites_to_trim_gaps_based = np.where(self.site_gappyness > gap_threshold)[0]
            col_char_freqs = self.column_character_frequencies()
            site_classification_types = np.array(
                [
                    determine_site_classification_type(col_char_freq)
                    for col_char_freq in col_char_freqs
                ]
            )
            sites_to_trim_classification_based = np.where(
                (
                    site_classification_types
                    != SiteClassificationType.parsimony_informative
                )
                | (site_classification_types != SiteClassificationType.constant)
            )[0]
            sites_to_trim = np.unique(
                np.concatenate(
                    (sites_to_trim_gaps_based, sites_to_trim_classification_based)
                )
            )[0]
        return sites_to_trim

    @staticmethod
    def from_bio_msa(alignment: MultipleSeqAlignment) -> "MSA":
        header_info = [
            {"id": rec.id, "name": rec.name, "description": rec.description}
            for rec in alignment
        ]
        seq_records = np.array([list(rec) for rec in alignment])
        return MSA(header_info, seq_records)

    def to_bio_msa(self) -> MultipleSeqAlignment:
        return self._to_bio_msa(self.sites_kept)

    def complement_to_bio_msa(self) -> MultipleSeqAlignment:
       return self._to_bio_msa(self.sites_trimmed)

    def _to_bio_msa(self, sites) -> MultipleSeqAlignment:
        return MultipleSeqAlignment(
            [
                SeqRecord(Seq("".join(rec)), **info)
                for rec, info in zip(sites.tolist(), self.header_info)
            ]
        )

    @property
    def stats(self) -> TrimmingStats:
        return TrimmingStats(self)

    def is_any_entry_sequence_only_gaps(
        self, gap_chars=DEFAULT_AA_GAP_CHARS
    ) -> tuple[bool, Union[str, None]]:
        # TODO: implement this w/ Jacob
        return False, None 
        # for entry, sequence in self._data.items():
        #     first_sequence_value = sequence[0]
        #     sequence_values_all_same = np.all(sequence == first_sequence_value)
        #     if (
        #         sequence_values_all_same
        #         and first_sequence_value.decode("utf-8") in gap_chars
        #     ):
        #         return True, entry

        # return False, None

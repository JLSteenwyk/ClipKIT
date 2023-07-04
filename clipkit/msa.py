from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from typing import Union

from .settings import DEFAULT_AA_GAP_CHARS


class MSA:
    def __init__(self, entries: list[str], starting_length: int) -> None:
        self.starting_length = starting_length
        self.entries = entries
        self._data = self._init_data()
        self._joined = None

    def _init_data(self) -> dict[str, np.array]:
        data = {}
        for entry in self.entries:
            data[entry] = np.zeros([self.starting_length], dtype=bytes)
        return data

    def _reset_joined(self) -> None:
        self._joined = None

    @property
    def length(self) -> int:
        return np.count_nonzero(next(iter(self._data.values())))

    @property
    def is_empty(self) -> bool:
        all_zeros = np.all(self._data[self.entries[0]] == b"")
        return all_zeros

    @property
    def string_sequences(self) -> dict[str, str]:
        if self._joined:
            return self._joined
        else:
            joined = {}
            for entry, sequence in self._data.items():
                joined[entry] = "".join(np.char.decode(sequence))
            self._joined = joined
        return self._joined

    @staticmethod
    def from_bio_msa(alignment: MultipleSeqAlignment) -> "MSA":
        entries = [entry.description for entry in alignment]
        length = alignment.get_alignment_length()
        return MSA(entries, length)

    def to_bio_msa(self) -> MultipleSeqAlignment:
        seqList = []
        for entry, joined in self.string_sequences.items():
            seqList.append(SeqRecord(Seq(str(joined)), id=str(entry), description=""))
        return MultipleSeqAlignment(seqList)

    def set_entry_sequence_at_position(
        self, entry: str, position: int, value: str
    ) -> None:
        self._data[entry][position] = value

        # since we've mutated the sequence data we need to release our cached _joined
        self._reset_joined()

    def is_any_entry_sequence_only_gaps(
        self, gap_chars=DEFAULT_AA_GAP_CHARS
    ) -> tuple[bool, Union[str, None]]:
        for entry, sequence in self._data.items():
            first_sequence_value = sequence[0]
            sequence_values_all_same = np.all(sequence == first_sequence_value)
            if (
                sequence_values_all_same
                and first_sequence_value.decode("utf-8") in gap_chars
            ):
                return True, entry

        return False, None

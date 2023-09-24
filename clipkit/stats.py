from dataclasses import dataclass
from typing import TYPE_CHECKING

from Bio.Align import MultipleSeqAlignment

if TYPE_CHECKING:
    from .msa import MSA


@dataclass
class TrimmingStats:
    msa: "MSA"

    @property
    def alignment_length(self) -> int:
        return self.msa.original_length

    @property
    def output_length(self) -> int:
        return self.msa.length

    @property
    def trimmed_length(self) -> int:
        return self.alignment_length - self.output_length

    @property
    def trimmed_percentage(self) -> float:
        return round((self.trimmed_length / self.alignment_length) * 100, 3)

    @property
    def summary(self) -> dict:
        return {
            "alignment_length": self.alignment_length,
            "output_length": self.output_length,
            "trimmed_length": self.trimmed_length,
            "trimmed_percentage": self.trimmed_percentage,
        }

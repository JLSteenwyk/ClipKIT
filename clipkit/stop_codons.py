from dataclasses import dataclass
from typing import Optional

from .modes import StopCodonMode


STOP_CODONS = frozenset({"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"})


@dataclass(frozen=True)
class StopCodonMaskingStats:
    mode: Optional[StopCodonMode] = None
    terminal_masked: int = 0
    internal_masked: int = 0

    @property
    def total_masked(self) -> int:
        return self.terminal_masked + self.internal_masked

    @property
    def summary(self) -> dict:
        return {
            "mode": self.mode.value if self.mode is not None else None,
            "terminal_masked": self.terminal_masked,
            "internal_masked": self.internal_masked,
            "total_masked": self.total_masked,
        }

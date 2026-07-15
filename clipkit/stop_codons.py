from dataclasses import dataclass
from typing import Optional, Union

from .exceptions import StopCodonValidationError
from .modes import StopCodonMode

STOP_CODONS = frozenset({"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"})


def normalize_stop_codon_mode(
    mode: Union[StopCodonMode, str, None],
) -> Optional[StopCodonMode]:
    if mode is None or isinstance(mode, StopCodonMode):
        return mode

    try:
        return StopCodonMode(mode)
    except ValueError as exc:
        choices = ", ".join(f"'{choice.value}'" for choice in StopCodonMode)
        raise StopCodonValidationError(
            f"remove_stop_codons must be one of: {choices}, or None."
        ) from exc


def validate_stop_codon_configuration(
    mode: Optional[StopCodonMode],
    *,
    codon: bool,
    sequence_type,
    alignment_length: int,
) -> None:
    if mode is None:
        return
    if not codon:
        raise StopCodonValidationError(
            "Stop codon masking requires codon-aware trimming "
            "(codon=True / --codon)."
        )
    if getattr(sequence_type, "value", sequence_type) != "nt":
        raise StopCodonValidationError(
            "Stop codon masking requires nucleotide input "
            "(sequence_type='nt' / --sequence_type nt)."
        )
    if alignment_length % 3 != 0:
        raise StopCodonValidationError(
            "Stop codon masking requires an alignment length divisible by 3."
        )


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

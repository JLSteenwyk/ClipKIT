"""IUPAC ambiguity definitions and helpers used during alignment analysis."""

from collections.abc import Iterable

NT_RESOLVED_STATES = frozenset("ACGTU")
AA_RESOLVED_STATES = frozenset("ACDEFGHIKLMNPQRSTVWYOU")

# Nucleotide ambiguity definitions use T internally. For RNA-only alignments,
# ``ambiguity_map`` substitutes U so fractional weights remain in the observed
# alphabet rather than introducing thymine into an RNA alignment.
_NT_AMBIGUITY_MAP = {
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
    "X": frozenset("ACGT"),
}

_AA_STANDARD_STATES = frozenset("ACDEFGHIKLMNPQRSTVWY")
_AA_AMBIGUITY_MAP = {
    "B": frozenset("DN"),
    "Z": frozenset("EQ"),
    "J": frozenset("IL"),
    "X": _AA_STANDARD_STATES,
}


def normalize_sequence_type(sequence_type: object) -> str:
    """Return ``aa`` or ``nt`` for enum and string sequence-type values."""
    value = getattr(sequence_type, "value", sequence_type)
    value = str(value).lower()
    if value not in {"aa", "nt"}:
        raise ValueError("sequence_type must be 'aa' or 'nt'")
    return value


def ambiguity_map(
    sequence_type: object, observed_states: Iterable[str] = ()
) -> dict[str, frozenset[str]]:
    """Return the relevant IUPAC ambiguity expansion map."""
    if normalize_sequence_type(sequence_type) == "aa":
        return _AA_AMBIGUITY_MAP

    observed_upper = {str(state).upper() for state in observed_states}
    if "U" not in observed_upper or "T" in observed_upper:
        return _NT_AMBIGUITY_MAP

    return {
        symbol: frozenset("U" if state == "T" else state for state in states)
        for symbol, states in _NT_AMBIGUITY_MAP.items()
    }


def ambiguity_symbols(sequence_type: object) -> frozenset[str]:
    """Return recognized ambiguity symbols for a sequence type."""
    if normalize_sequence_type(sequence_type) == "aa":
        return frozenset(_AA_AMBIGUITY_MAP)
    return frozenset(_NT_AMBIGUITY_MAP)


def resolved_symbols(sequence_type: object) -> frozenset[str]:
    """Return unambiguous biological states for a sequence type."""
    if normalize_sequence_type(sequence_type) == "aa":
        return AA_RESOLVED_STATES
    return NT_RESOLVED_STATES


def is_nucleotide_alphabet(characters: Iterable[str]) -> bool:
    """Whether all observed characters belong to the IUPAC nucleotide alphabet."""
    nucleotide_symbols = NT_RESOLVED_STATES | ambiguity_symbols("nt")
    ignored_symbols = frozenset("-?*.")
    observed = {str(character).upper() for character in characters}
    return observed.issubset(nucleotide_symbols | ignored_symbols)

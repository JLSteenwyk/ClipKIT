import logging

from .helpers import SeqType
from .settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .msa import MSA


logger = logging.getLogger(__name__)


def warn_if_all_sites_were_trimmed(msa: "MSA") -> None:
    if msa.is_empty:
        logger.warning(
            "WARNING: All sites trimmed from alignment. Please use different parameters."
        )


def warn_if_entry_contains_only_gaps(msa: "MSA") -> None:
    should_warn, entry = msa.is_any_entry_sequence_only_gaps()
    if should_warn:
        logger.warning(f"WARNING: header id '{entry}' contains only gaps")

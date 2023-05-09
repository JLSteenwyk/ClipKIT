import logging

from .helpers import SeqType
from .settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS

logger = logging.getLogger(__name__)


def warn_if_all_sites_were_trimmed(keepMSA):
    if keepMSA.is_empty:
        logger.warning(
            "WARNING: All sites trimmed from alignment. Please use different parameters."
        )


def warn_if_entry_contains_only_gaps(keepMSA, sequence_type):
    if sequence_type == SeqType.aa:
        gap_chars = DEFAULT_AA_GAP_CHARS
    else:
        gap_chars = DEFAULT_NT_GAP_CHARS
    should_warn, entry = keepMSA.is_any_entry_sequence_only_gaps(gap_chars)
    if should_warn:
        logger.warning(f"WARNING: header id '{entry}' contains only gaps")

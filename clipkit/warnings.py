import logging

logger = logging.getLogger(__name__)


def warn_if_all_sites_were_trimmed(keepD):
    if not len(next(iter(keepD.values()))):
        logger.warning(
            "WARNING: All sites trimmed from alignment. Please use different parameters."
        )


def warn_if_entry_contains_only_gaps(keepD):
    for entry, sequence in keepD.items():
        chars_in_sequence = set(sequence)
        if len(chars_in_sequence) == 1 and "-" in chars_in_sequence:
            logger.warning(f"WARNING: header id '{entry}' contains only gaps")
            break

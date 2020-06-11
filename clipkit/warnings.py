import logging

logger = logging.getLogger(__name__)


def check_if_all_sites_were_trimmed(keepD):
    """
    Determine if every site in an alignment was trimmed
    """
    # check if resulting alingment length is 0
    if not len(next(iter(keepD.values()))):
        logger.warning(
            "WARNING: All sites trimmed from alignment. Please use different parameters."
        )


def check_if_entry_contains_only_gaps(keepD):
    """
    Determine if any output sequence contains only gaps
    """
    for entry, sequence in keepD.items():
        chars_in_sequence = set(sequence)
        if len(chars_in_sequence) == 1 and "-" in chars_in_sequence:
            logger.warning(f"WARNING: header id '{entry}' contains only gaps")
            break

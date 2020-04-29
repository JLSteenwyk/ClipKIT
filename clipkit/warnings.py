import logging

logger = logging.getLogger(__name__)

## Function to determine if every site in an alignment was trimmed
def checking_if_all_sites_were_trimmed(keepD):
    # check if resulting alingment length is 0
    if not len(next(iter(keepD.values()))):
        logger.warning("WARNING: All sites trimmed from alignment. Please use different parameters.")

## Function to determine if any output sequence contains only gaps
def checking_if_entry_contains_only_gaps(keepD):
    # checking if any sequence entry contains only gaps
    for entry, sequence in keepD.items():
        chars_in_sequence = set(sequence)
        if len(chars_in_sequence) == 1 and '-' in chars_in_sequence:
            logger.warning(f"WARNING: header id '{entry}' contains only gaps")
            break


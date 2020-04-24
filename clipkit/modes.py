import logging

from enum import Enum

logger = logging.getLogger("clipkit.clipkit")
####################################################################
### Mode Functions 	                                             ###
### This block of code contains all functions that populate      ###
### dictionaries with sites to trim and keep according to user   ###
### set mode                                           			 ###
####################################################################


class TrimmingMode(Enum):
    kpi_gappy = "kpi-gappy"
    kpi = "kpi"
    gappy = "gappy"

# TODO: write unit test
def shouldKeep(mode, parsimony_informative, gappyness, gaps):
    print(f"mode: {mode}, parsimony_informative: {parsimony_informative}, gappyness: {gappyness}, gaps: {gaps}")
    if mode == TrimmingMode.kpi_gappy:
        return gappyness <= gaps and parsimony_informative
    elif mode == TrimmingMode.gappy:
        return gappyness <= gaps
    elif mode == TrimmingMode.kpi:
        return parsimony_informative

def trim(
    gappyness,
    parsimony_informative,
    keepD,
    trimD,
    alignment_position: int,
    gaps,
    alignment,
    mode: TrimmingMode,
    use_log: bool,
):
    """
    Trims according to the mode kpi-gappy wherein only parismony informative
    sites are kept. Additionally, sites that are sufficiently gappy are removed

    Arguments
    ---------
    argv: gappyness
        a float value for what fraction of characters are gaps
    argv: parsimony_informative
        boolean for whether or not a site is parsimony informative
    argv: keepD
        dictionary for sites that will be kept in the resulting alignment
    argv: trimD
        dictionary for sites that will be removed from the final alignment
    argv: i
        i is the position in the alignment that the loop is currently in
    argv: gaps
        gaps threshold to determine if a position is too gappy or not
    argv: alignment
        biopython multiple sequence alignment object
    """
    # save to keepD
    if shouldKeep(mode, parsimony_informative, gappyness, gaps):
        for entry in alignment:
            keepD[entry.id][alignment_position] = (entry.seq._data[alignment_position])
        if use_log:
            logger.debug(f"{str(alignment_position + 1)} keep PI {gappyness}")
    else:
        # save to trimD
        for entry in alignment:
            trimD[entry.id][alignment_position] = (entry.seq._data[alignment_position])
        if use_log:
            logger.debug(f"{str(alignment_position + 1)} trim nPI {gappyness}")

    return keepD, trimD

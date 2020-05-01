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
    kpic = "kpic"
    kpic_gappy = "kpic-gappy"

# Function to determine if a site should be kept or not
def shouldKeep(mode, parsimony_informative, constant_site, gappyness, gaps):
    # print(f"mode: {mode}, parsimony_informative: {parsimony_informative}, gappyness: {gappyness}, gaps: {gaps}")
    if mode == TrimmingMode.kpi_gappy:
        return gappyness <= gaps and parsimony_informative
    elif mode == TrimmingMode.gappy:
        return gappyness <= gaps
    elif mode == TrimmingMode.kpi:
        return parsimony_informative
    elif mode == TrimmingMode.kpic:
        return parsimony_informative or constant_site
    elif mode == TrimmingMode.kpic_gappy:
        return gappyness <= gaps and (parsimony_informative or constant_site)

def trim(
    gappyness: float,
    parsimony_informative: bool,
    constant_site: bool,
    keepD: dict,
    trimD: dict,
    alignment_position: int,
    gaps: float,
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
    if shouldKeep(mode, parsimony_informative, constant_site, gappyness, gaps):
        for entry in alignment:
            keepD[entry.id][alignment_position] = (entry.seq._data[alignment_position])
        if use_log:
            if constant_site:
                logger.debug(f"{str(alignment_position + 1)} keep Const {gappyness}")
            elif parsimony_informative:
                logger.debug(f"{str(alignment_position + 1)} keep PI {gappyness}")
            else:
                logger.debug(f"{str(alignment_position + 1)} keep nConst,nPI {gappyness}")
    # save to trimD
    else:
        for entry in alignment:
            trimD[entry.id][alignment_position] = (entry.seq._data[alignment_position])
        if use_log:
            if constant_site:
                logger.debug(f"{str(alignment_position + 1)} trim Const {gappyness}")
            elif parsimony_informative:
                logger.debug(f"{str(alignment_position + 1)} trim PI {gappyness}")
            else:
                logger.debug(f"{str(alignment_position + 1)} trim nConst,nPI {gappyness}")

    return keepD, trimD

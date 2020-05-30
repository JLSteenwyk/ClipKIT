import logging

from enum import Enum

logger = logging.getLogger("clipkit.clipkit")


class TrimmingMode(Enum):
    kpi_gappy = "kpi-gappy"
    kpi = "kpi"
    gappy = "gappy"
    kpic = "kpic"
    kpic_gappy = "kpic-gappy"
    medium = "medium"
    medium_gappy = "medium-gappy"
    heavy = "heavy"
    heavy_gappy = "heavy-gappy"


def shouldKeep(
    mode: TrimmingMode,
    parsimony_informative: bool,
    constant_site: bool,
    gappyness: float,
    gaps: float,
):
    """
    Determine if a site should be kept or not
    """
    if mode in (TrimmingMode.kpi_gappy, TrimmingMode.heavy_gappy):
        return gappyness <= gaps and parsimony_informative
    elif mode == TrimmingMode.gappy:
        return gappyness <= gaps
    elif mode in (TrimmingMode.kpi, TrimmingMode.heavy):
        return parsimony_informative
    elif mode in (TrimmingMode.kpic, TrimmingMode.medium):
        return parsimony_informative or constant_site
    elif mode in (TrimmingMode.kpic_gappy, TrimmingMode.medium_gappy):
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
            keepD[entry.id][alignment_position] = entry.seq._data[alignment_position]
        if use_log:
            if constant_site:
                logger.debug(f"{str(alignment_position + 1)} keep Const {gappyness}")
            elif parsimony_informative:
                logger.debug(f"{str(alignment_position + 1)} keep PI {gappyness}")
            else:
                logger.debug(
                    f"{str(alignment_position + 1)} keep nConst,nPI {gappyness}"
                )
    # save to trimD
    else:
        for entry in alignment:
            trimD[entry.id][alignment_position] = entry.seq._data[alignment_position]
        if use_log:
            if constant_site:
                logger.debug(f"{str(alignment_position + 1)} trim Const {gappyness}")
            elif parsimony_informative:
                logger.debug(f"{str(alignment_position + 1)} trim PI {gappyness}")
            else:
                logger.debug(
                    f"{str(alignment_position + 1)} trim nConst,nPI {gappyness}"
                )

    return keepD, trimD

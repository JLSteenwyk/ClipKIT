import logging

from enum import Enum

logger = logging.getLogger('clipkit.clipkit')
####################################################################
### Mode Functions 	                                             ###
### This block of code contains all functions that populate      ###
### dictionaries with sites to trim and keep according to user   ###
### set mode                                           			 ###
####################################################################


class TrimmingMode(Enum):
    kpi_gappy = 'kpi-gappy'
    kpi = 'kpi'
    gappy = 'gappy'

def kpi_gappy_mode(
    gappyness, parsimony_informative, 
    keepD, trimD, i, gaps, alignment
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

    # if gappyness is lower than gaps threshold and the site is
    # parismony informative, save to keepD. Otherwise, save to trimD.
    # All the while, save information to logD
    if gappyness <= gaps and parsimony_informative:
        # save to keepD
        for entry in alignment:
            keepD[entry.id].append(entry.seq._data[i])
        # save to logL - structure is site, parsimony informative (PI) or not (nPI)
        # gappyness, kept or trimmed
        temp = []
        temp.append(str(i+1))
        temp.append('keep')
        temp.append("PI")
        temp.append(str(gappyness))
        logger.info(temp)
    else:
        # save to trimD 
        for entry in alignment:
            trimD[entry.id].append(entry.seq._data[i])
        # save to logL - structure is site, parsimony informative (PI) or not (nPI)
        # gappyness, kept or trimmed
        temp = []
        temp.append(str(i+1))
        temp.append('trim')
        temp.append("nPI")
        temp.append(str(gappyness))
        logger.info(temp)
    
    return keepD, trimD

def gappy_mode(
    gappyness, parsimony_informative, 
    keepD, trimD, i, gaps, alignment
    ):
    """
    Trims according to the mode gappy wherein sites that are sufficiently gappy 
    are removed

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

    # if gappyness is lower than gaps threshold and the site is
    # parismony informative, save to keepD. Otherwise, save to trimD.
    # All the while, save information to logD
    if gappyness <= gaps:
        # save to keepD
        for entry in alignment:
            keepD[entry.id].append(entry.seq._data[i])
        # save to logL - structure is site, parsimony informative (PI) or not (nPI)
        # gappyness, kept or trimmed
        if parsimony_informative:
            temp = []
            temp.append(str(i+1))
            temp.append('keep')
            temp.append("PI")
            temp.append(str(gappyness))
            logger.info(temp)
        else: 
            temp = []
            temp.append(str(i+1))
            temp.append('keep')
            temp.append("nPI")
            temp.append(str(gappyness))
            logger.info(temp)
    else:
        # save to trimD 
        for entry in alignment:
            trimD[entry.id].append(entry.seq._data[i])
        # save to logL - structure is site, parsimony informative (PI) or not (nPI)
        # gappyness, kept or trimmed
        if parsimony_informative:
            temp = []
            temp.append(str(i+1))
            temp.append('trim')
            temp.append("PI")
            temp.append(str(gappyness))
            logger.info(temp)
        else: 
            temp = []
            temp.append(str(i+1))
            temp.append('trim')
            temp.append("nPI")
            temp.append(str(gappyness))
            logger.info(temp)
    
    return keepD, trimD

def kpi_mode(
    gappyness, parsimony_informative, 
    keepD, trimD, i, gaps, alignment
    ):
    """
    Trims according to the mode gappy wherein sites that are sufficiently gappy 
    are removed

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
    # if site is parismony informative, save to keepD. 
    # Otherwise, save to trimD.
    # All the while, save information to logD
    if parsimony_informative:
        # save to keepD
        for entry in alignment:
            keepD[entry.id].append(entry.seq._data[i])

        temp = []
        temp.append(str(i+1))
        temp.append('keep')
        temp.append("PI")
        temp.append(str(gappyness))
        logger.info(temp)
    else:
        # save to trimD 
        for entry in alignment:
            trimD[entry.id].append(entry.seq._data[i])

        temp = []
        temp.append(str(i+1))
        temp.append('trim')
        temp.append("nPI")
        temp.append(str(gappyness))
        logger.info(temp)

    return keepD, trimD
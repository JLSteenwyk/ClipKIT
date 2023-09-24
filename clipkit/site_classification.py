from enum import Enum


class SiteClassificationType(Enum):
    parsimony_informative = "parsimony-informative"
    constant = "constant"
    singleton = "singleton"
    other = "other"


def determine_site_classification_type(
    character_counts: dict,
) -> SiteClassificationType:
    """
    Determines if a site is parsimony informative or constant.
    A site is parsimony-informative if it contains at least two types of nucleotides
    (or amino acids), and at least two of them occur with a minimum frequency of two.
    https://www.megasoftware.net/web_help_7/rh_parsimony_informative_site.htm

    A site is constant if it contains only one character and that character occurs
    at least twice. https://www.megasoftware.net/web_help_7/rh_constant_site.htm

    A singleton is a site that contains at least two types of characters with, at most,
    one occuring multiple times. https://www.megasoftware.net/web_help_7/rh_singleton_sites.htm
    """
    parsimony_informative_threshold = 2
    counts_gte_threshold = 0

    for count in character_counts.values():
        if count >= 2:
            counts_gte_threshold += 1
        if counts_gte_threshold >= parsimony_informative_threshold:
            return SiteClassificationType.parsimony_informative

    if counts_gte_threshold == 1 and len(character_counts) == 1:
        return SiteClassificationType.constant
    elif counts_gte_threshold == 1 and len(character_counts) > 1:
        return SiteClassificationType.singleton

    return SiteClassificationType.other

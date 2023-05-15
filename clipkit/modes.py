from .logger import log_file_logger
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .helpers import SiteClassificationType


class TrimmingMode(Enum):
    gappy = "gappy"
    smart_gap = "smart-gap"
    kpi = "kpi"
    kpi_gappy = "kpi-gappy"
    kpi_smart_gap = "kpi-smart-gap"
    kpic = "kpic"
    kpic_gappy = "kpic-gappy"
    kpic_smart_gap = "kpic-smart-gap"


def shouldKeep(
    mode: TrimmingMode,
    site_classification_type: "SiteClassificationType",
    gappyness: float,
    gaps: float,
):
    """
    Determine if a site should be kept or not
    """
    if mode == TrimmingMode.kpi_gappy:
        return gappyness <= gaps and site_classification_type.parsimony_informative
    elif mode == TrimmingMode.gappy:
        return gappyness <= gaps
    elif mode == TrimmingMode.kpi:
        return site_classification_type.parsimony_informative
    elif mode == TrimmingMode.kpic:
        return (
            site_classification_type.parsimony_informative
            or site_classification_type.constant_site
        )
    elif mode == TrimmingMode.kpic_gappy:
        return gappyness <= gaps and (
            site_classification_type.parsimony_informative
            or site_classification_type.constant_site
        )
    elif mode == TrimmingMode.smart_gap:
        return round(gappyness, 4) < gaps
    elif mode == TrimmingMode.kpic_smart_gap:
        return round(gappyness, 4) < gaps and (
            site_classification_type.parsimony_informative
            or site_classification_type.constant_site
        )
    elif mode == TrimmingMode.kpi_smart_gap:
        return (
            round(gappyness, 4) < gaps
            and site_classification_type.parsimony_informative
        )


def trim(
    gappyness: float,
    site_classification_type: "SiteClassificationType",
    keepMSA: dict,
    trimMSA: dict,
    alignment_position: int,
    gaps: float,
    alignment,
    mode: TrimmingMode,
    use_log: bool,
):
    if shouldKeep(mode, site_classification_type, gappyness, gaps):
        for entry in alignment:
            new_value = entry.seq._data[alignment_position : alignment_position + 1]
            keepMSA.set_entry_sequence_at_position(
                entry.description, alignment_position, new_value
            )
        if use_log:
            log_file_logger.debug(
                f"{str(alignment_position + 1)} keep {site_classification_type.value} {gappyness}"
            )
    elif trimMSA is not None:
        for entry in alignment:
            new_value = entry.seq._data[alignment_position : alignment_position + 1]
            trimMSA.set_entry_sequence_at_position(
                entry.description, alignment_position, new_value
            )
        if use_log:
            log_file_logger.debug(
                f"{str(alignment_position + 1)} trim {site_classification_type.value} {gappyness}"
            )

    return keepMSA, trimMSA

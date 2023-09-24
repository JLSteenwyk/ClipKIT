from enum import Enum
from typing import TYPE_CHECKING
from .logger import log_file_logger

if TYPE_CHECKING:
    from Bio.Align import MultipleSeqAlignment
    from .msa import MSA


class TrimmingMode(Enum):
    gappy = "gappy"
    smart_gap = "smart-gap"
    kpi = "kpi"
    kpi_gappy = "kpi-gappy"
    kpi_smart_gap = "kpi-smart-gap"
    kpic = "kpic"
    kpic_gappy = "kpic-gappy"
    kpic_smart_gap = "kpic-smart-gap"

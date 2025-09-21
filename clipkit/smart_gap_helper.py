from collections import Counter
from functools import lru_cache

from Bio.Align import MultipleSeqAlignment
import numpy as np

from .logger import logger


@lru_cache(maxsize=128)
def _cached_gap_slope(gap_tuple, alignment_length):
    """Cached calculation of gap-to-gap slope for a specific gap distribution"""
    gaps_arr = np.array(gap_tuple)
    return gap_to_gap_slope(gaps_arr, alignment_length)


def smart_gap_threshold_determination(
    alignment: MultipleSeqAlignment, gap_chars: list
) -> float:
    alignment_length = alignment.get_alignment_length()

    # get distribution of gaps rounded to the fourth decimal place
    gaps_dist = get_gaps_distribution_optimized(alignment, gap_chars)

    # count freq of gaps and convert to sorted np array
    gaps_arr = count_and_sort_gaps_optimized(gaps_dist)

    # calculate gap-to-gap slope
    slopes = gap_to_gap_slope_vectorized(gaps_arr, alignment_length)

    # find the greatest difference in slopes and set gaps to y1
    return greatest_diff_in_slopes(slopes, gaps_arr)


def greatest_diff_in_slopes(slopes: list[float], gaps_arr: np.array) -> float:
    diffs = []
    # if there is only one slope, use that value to determine
    # the threshold. Otherwise, calculate the greatest difference
    # in slopes
    if len(slopes) > 1:
        # Vectorized difference calculation
        slopes_arr = np.array(slopes)
        diffs = np.abs(np.diff(slopes_arr))
        max_diff_idx = np.argmax(diffs)
        return gaps_arr[max_diff_idx][0]
    elif len(slopes) == 0:
        return 1
    else:
        return gaps_arr[0][0]


def gap_to_gap_slope_vectorized(gaps_arr: np.array, alignment_length: int) -> list[float]:
    """Vectorized calculation of gap-to-gap slopes"""
    if len(gaps_arr) < 2:
        return []

    # Convert to numpy array for vectorized operations
    gaps_values = gaps_arr[:, 0].astype(float)
    gaps_counts = gaps_arr[:, 1].astype(float)

    # Normalize counts by alignment length
    normalized_counts = gaps_counts / alignment_length

    # Calculate cumulative sums
    cumsum = np.cumsum(normalized_counts)

    # Calculate slopes between consecutive points
    slopes = []
    for i in range(1, len(gaps_arr)):
        if gaps_values[i] != gaps_values[i-1]:  # Avoid division by zero
            slope = abs((cumsum[i] - cumsum[i-1]) / (gaps_values[i] - gaps_values[i-1]))
            slopes.append(slope)

    # Only use first half of slopes
    return slopes[: (len(slopes) // 2)]


def gap_to_gap_slope(gaps_arr: np.array, alignment_length: int) -> list[float]:
    """Original implementation kept for compatibility"""
    sum_sites_current = gaps_arr[0][1] / alignment_length
    sum_sites_previous = 0
    slopes = []
    for previous, current in zip(gaps_arr, gaps_arr[1:]):
        sum_sites_previous += previous[1] / alignment_length
        sum_sites_current += current[1] / alignment_length
        slopes.append(
            abs((sum_sites_current - sum_sites_previous) / (current[0] - previous[0]))
        )
    # only use first half of slopes
    return slopes[: (len(slopes) // 2)]


def get_gaps_distribution_optimized(
    alignment: MultipleSeqAlignment, gap_chars: list
) -> np.ndarray:
    """Optimized gap distribution calculation using numpy broadcasting"""
    # Convert alignment to numpy array once
    msa_array = np.array([list(rec) for rec in alignment], dtype='U1')

    # Use broadcasting for gap checking
    gap_mask = np.zeros(msa_array.shape, dtype=bool)
    for gap_char in gap_chars:
        gap_mask |= (msa_array == gap_char)

    # Calculate mean across sequences (axis 0)
    gaps_dist = gap_mask.mean(axis=0)

    return np.round(gaps_dist, decimals=4)


def get_gaps_distribution(
    alignment: MultipleSeqAlignment, gap_chars: list
) -> list[float]:
    """Original implementation kept for compatibility"""
    msa_array = np.array([list(rec) for rec in alignment])
    gaps_dist = (np.isin(msa_array, gap_chars)).mean(axis=0)

    return np.round(gaps_dist, decimals=4).tolist()


def count_and_sort_gaps_optimized(gaps_dist: np.ndarray) -> np.ndarray:
    """Optimized counting and sorting using numpy operations"""
    # Use numpy's unique for counting
    unique_gaps, counts = np.unique(gaps_dist, return_counts=True)

    # Stack into 2D array
    gaps_arr = np.column_stack((unique_gaps, counts))

    # Sort by gap value (descending)
    sorted_indices = np.argsort(-gaps_arr[:, 0])

    return gaps_arr[sorted_indices]


def count_and_sort_gaps(gaps_dist: list) -> list:
    """Original implementation kept for compatibility"""
    gaps_count = dict(Counter(gaps_dist))
    gaps_arr = np.array(list(gaps_count.items()))
    return gaps_arr[np.argsort(-gaps_arr[:, 0])].tolist()
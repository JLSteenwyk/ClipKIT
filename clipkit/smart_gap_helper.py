from collections import Counter
from typing import Optional

from Bio.Align import MultipleSeqAlignment
import numpy as np


def smart_gap_threshold_determination(
    alignment: MultipleSeqAlignment,
    gap_chars: list,
    seq_records: Optional[np.ndarray] = None,
    gaps_dist: Optional[np.ndarray] = None,
) -> float:
    alignment_length = (
        seq_records.shape[1]
        if seq_records is not None
        else alignment.get_alignment_length()
    )

    # get distribution of gaps rounded to the fourth decimal place
    if gaps_dist is None:
        gaps_dist = get_gaps_distribution_optimized(
            alignment, gap_chars, seq_records=seq_records
        )

    # count freq of gaps and convert to sorted np array
    gaps_arr = count_and_sort_gaps_optimized(gaps_dist)

    # calculate gap-to-gap slope
    slopes = gap_to_gap_slope_vectorized(gaps_arr, alignment_length)

    # find the greatest difference in slopes and set gaps to y1
    return greatest_diff_in_slopes(slopes, gaps_arr)


def greatest_diff_in_slopes(slopes: list[float], gaps_arr: np.array) -> float:
    if len(slopes) == 0:
        return 1
    if len(slopes) == 1:
        return gaps_arr[0][0]
    max_diff_idx = np.argmax(np.abs(np.diff(np.array(slopes))))
    return gaps_arr[max_diff_idx][0]


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

    distinct = gaps_values[1:] != gaps_values[:-1]
    gap_changes = gaps_values[1:][distinct] - gaps_values[:-1][distinct]
    count_changes = cumsum[1:][distinct] - cumsum[:-1][distinct]
    slopes = np.abs(count_changes / gap_changes)
    return list(slopes[: len(slopes) // 2])


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
    alignment: MultipleSeqAlignment,
    gap_chars: list,
    seq_records: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Optimized gap distribution calculation using numpy broadcasting"""
    # Reuse precomputed sequence matrix when available.
    msa_array = (
        seq_records
        if seq_records is not None
        else np.array([list(rec) for rec in alignment], dtype="U1")
    )

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

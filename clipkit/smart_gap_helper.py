from collections import Counter

from Bio.Align import MultipleSeqAlignment
import numpy as np
from tqdm import tqdm

from .helpers import report_column_features
from .logger import logger


def smart_gap_threshold_determination(
    alignment: MultipleSeqAlignment, gap_chars: list, quiet: bool
) -> float:
    alignment_length = alignment.get_alignment_length()

    # get distribution of gaps rounded to the fourth decimal place
    gaps_dist = get_gaps_distribution(alignment, alignment_length, gap_chars, quiet)

    # count freq of gaps and convert to sorted np array
    gaps_arr = count_and_sort_gaps(gaps_dist)

    # calculate gap-to-gap slope
    slopes = gap_to_gap_slope(gaps_arr, alignment_length)

    # find the greatest difference in slopes and set gaps to y1
    return greatest_diff_in_slopes(slopes, gaps_arr)


def greatest_diff_in_slopes(slopes: list[float], gaps_arr: np.array) -> float:
    diffs = []
    # if there is only one slope, use that value to determine
    # the threshold. Otherwise, calculate the greatest difference
    # in slopes
    if len(slopes) > 1:
        for val0, val1 in zip(slopes, slopes[1:]):
            diff0 = abs(val0 - val1)
            diffs.append(diff0)
    elif len(slopes) == 0:
        return 1
    else:
        diffs = slopes
    return gaps_arr[diffs.index(max(diffs))][0]


def gap_to_gap_slope(gaps_arr: np.array, alignment_length: int) -> list[float]:
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


def get_gaps_distribution(
    alignment: MultipleSeqAlignment, alignment_length: int, gap_chars: list, quiet: bool
) -> list[float]:
    gaps_dist = []
    for i in tqdm(range(0, alignment_length), disable=quiet, postfix="smart-gap"):
        _, gappyness = report_column_features(alignment, i, gap_chars)
        gappyness = round(gappyness, 4)
        gaps_dist.append(gappyness)

    return gaps_dist


def count_and_sort_gaps(gaps_dist: list) -> list:
    gaps_count = dict(Counter(gaps_dist))
    gaps_arr = np.array(list(gaps_count.items()))
    return gaps_arr[np.argsort(-gaps_arr[:, 0])].tolist()

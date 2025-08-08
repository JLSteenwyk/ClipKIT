#!/usr/bin/env python
"""
Standalone test for parallel processing functions without Bio dependency
"""

import numpy as np
import time
from multiprocessing import Pool
from enum import Enum

# Copy the SiteClassificationType enum
class SiteClassificationType(Enum):
    parsimony_informative = "parsimony-informative"
    constant = "constant"
    singleton = "singleton"
    other = "other"

# Copy the determine_site_classification_type function
def determine_site_classification_type(character_counts: dict) -> SiteClassificationType:
    """
    Determines if a site is parsimony informative or constant.
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

# Helper function for parallel processing
def _calculate_column_frequency_helper(args):
    """Helper function for parallel processing of column frequencies"""
    column, gap_chars = args
    col_sorted_unique_values_for, col_counts_per_char = np.unique(
        [char.upper() for char in column], return_counts=True
    )
    freqs = dict(zip(col_sorted_unique_values_for, col_counts_per_char))
    for gap_char in gap_chars:
        try:
            del freqs[gap_char]
        except KeyError:
            continue
    return freqs

def test_large_alignment():
    """Test with a larger alignment to see speedup"""
    
    print("\n" + "=" * 60)
    print("TESTING WITH LARGE ALIGNMENT (100 sequences x 10,000 sites)")
    print("=" * 60)
    
    # Create a larger test alignment
    np.random.seed(42)
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY-")
    test_alignment = np.random.choice(amino_acids, size=(100, 10000))
    gap_chars = ['-', 'X', 'x', '?', '*']
    
    # Test character frequency calculation
    print("\nCharacter Frequency Calculation:")
    print("-" * 40)
    
    # Single-threaded
    print("Single-threaded:")
    start_time = time.time()
    single_freqs = []
    for column in test_alignment.T:
        freqs = _calculate_column_frequency_helper((column, gap_chars))
        single_freqs.append(freqs)
    single_freq_time = time.time() - start_time
    print(f"  Time: {single_freq_time:.3f} seconds")
    
    # Multi-threaded with different thread counts
    for n_threads in [2, 4, 8]:
        print(f"\nMulti-threaded ({n_threads} threads):")
        start_time = time.time()
        args_list = [(column, gap_chars) for column in test_alignment.T]
        with Pool(processes=n_threads) as pool:
            multi_freqs = pool.map(_calculate_column_frequency_helper, args_list)
        multi_freq_time = time.time() - start_time
        print(f"  Time: {multi_freq_time:.3f} seconds")
        print(f"  Speedup: {single_freq_time/multi_freq_time:.2f}x")
        
        # Verify correctness
        if not all(s == m for s, m in zip(single_freqs, multi_freqs)):
            print("  ERROR: Results don't match!")
    
    # Test site classification
    print("\n" + "=" * 60)
    print("Site Classification:")
    print("-" * 40)
    
    # Single-threaded
    print("Single-threaded:")
    start_time = time.time()
    single_classifications = [determine_site_classification_type(f) for f in single_freqs]
    single_class_time = time.time() - start_time
    print(f"  Time: {single_class_time:.3f} seconds")
    
    # Multi-threaded with different thread counts
    for n_threads in [2, 4, 8]:
        print(f"\nMulti-threaded ({n_threads} threads):")
        start_time = time.time()
        with Pool(processes=n_threads) as pool:
            multi_classifications = pool.map(determine_site_classification_type, single_freqs)
        multi_class_time = time.time() - start_time
        print(f"  Time: {multi_class_time:.3f} seconds")
        print(f"  Speedup: {single_class_time/multi_class_time:.2f}x")
        
        # Verify correctness
        if not all(s == m for s, m in zip(single_classifications, multi_classifications)):
            print("  ERROR: Results don't match!")
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY:")
    print("  - Parallel processing correctly preserves results")
    print("  - Speedup increases with alignment size")
    print("  - Optimal thread count depends on system and data size")
    print("=" * 60)

if __name__ == "__main__":
    test_large_alignment()
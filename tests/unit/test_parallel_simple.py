#!/usr/bin/env python
"""
Simple test for parallel processing functions
"""

import numpy as np
import time
from multiprocessing import Pool

# Copy the helper function from msa.py
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

def test_parallel_frequency_calculation():
    """Test that parallel processing of character frequencies works"""
    
    # Create a test alignment (10 sequences, 1000 sites)
    np.random.seed(42)
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY-")
    test_alignment = np.random.choice(amino_acids, size=(10, 1000))
    gap_chars = ['-', 'X', 'x', '?', '*']
    
    print("Testing parallel character frequency calculation...")
    print("-" * 50)
    
    # Single-threaded version
    print("\n1. Single-threaded calculation:")
    start_time = time.time()
    single_results = []
    for column in test_alignment.T:
        freqs = _calculate_column_frequency_helper((column, gap_chars))
        single_results.append(freqs)
    single_time = time.time() - start_time
    print(f"   Completed in {single_time:.3f} seconds")
    
    # Multi-threaded version
    print("\n2. Multi-threaded calculation (4 threads):")
    start_time = time.time()
    args_list = [(column, gap_chars) for column in test_alignment.T]
    with Pool(processes=4) as pool:
        multi_results = pool.map(_calculate_column_frequency_helper, args_list)
    multi_time = time.time() - start_time
    print(f"   Completed in {multi_time:.3f} seconds")
    
    # Compare results
    print("\n3. Verifying results:")
    results_match = all(s == m for s, m in zip(single_results, multi_results))
    if results_match:
        print("   ✓ Results are identical - parallel processing preserves correctness")
    else:
        print("   ✗ Results differ - there may be an issue")
        assert results_match
    
    # Report speedup
    print("\n" + "-" * 50)
    if multi_time < single_time:
        speedup = single_time / multi_time
        print(f"Speedup with 4 threads: {speedup:.2f}x")
    else:
        print("Note: Parallel version was slower (overhead for small dataset)")
    
def test_site_classification():
    """Test site classification (from site_classification.py)"""
    from clipkit.site_classification import determine_site_classification_type, SiteClassificationType
    
    print("\n\nTesting parallel site classification...")
    print("-" * 50)
    
    # Generate test character frequencies
    test_freqs = [
        {'A': 5, 'T': 3},  # Parsimony informative
        {'A': 10},         # Constant
        {'A': 9, 'T': 1},  # Singleton
        {'A': 1, 'T': 1, 'G': 1}  # Other
    ] * 250  # Create 1000 test sites
    
    # Single-threaded
    print("\n1. Single-threaded classification:")
    start_time = time.time()
    single_results = [determine_site_classification_type(f) for f in test_freqs]
    single_time = time.time() - start_time
    print(f"   Completed in {single_time:.3f} seconds")
    
    # Multi-threaded
    print("\n2. Multi-threaded classification (4 threads):")
    start_time = time.time()
    with Pool(processes=4) as pool:
        multi_results = pool.map(determine_site_classification_type, test_freqs)
    multi_time = time.time() - start_time
    print(f"   Completed in {multi_time:.3f} seconds")
    
    # Verify
    print("\n3. Verifying results:")
    results_match = all(s == m for s, m in zip(single_results, multi_results))
    if results_match:
        print("   ✓ Results are identical")
    else:
        print("   ✗ Results differ")
        assert results_match
    
    print("-" * 50)
    if multi_time < single_time:
        speedup = single_time / multi_time
        print(f"Speedup with 4 threads: {speedup:.2f}x")
    
if __name__ == "__main__":
    print("=" * 50)
    print("CLIPKIT PARALLEL PROCESSING TEST")
    print("=" * 50)
    
    test_parallel_frequency_calculation()
    test_site_classification()
    
    print("\n" + "=" * 50)
    print("All tests passed!")
    print("=" * 50)

#!/usr/bin/env python
"""
Test script for parallel processing implementation in ClipKIT
"""

import time
import sys
from clipkit import clipkit as ck
from clipkit.files import FileFormat
from clipkit.modes import TrimmingMode
from clipkit.helpers import SeqType

def test_parallel_processing():
    """Test that parallel processing works correctly"""
    
    # Test with a sample alignment file
    test_file = "tests/integration/samples/EOG092C0CZK_aa_aln.fasta"
    
    print("Testing ClipKIT with parallel processing...")
    print("-" * 50)
    
    # Test with single thread
    print("\n1. Testing with single thread (--threads 1):")
    start_time = time.time()
    
    try:
        ck.execute(
            input_file=test_file,
            input_file_format=None,
            output_file="test_output_single.fa",
            output_file_format=FileFormat.fasta,
            sequence_type=None,
            gaps=0.9,
            gap_characters=None,
            complement=False,
            codon=False,
            ends_only=False,
            mode=TrimmingMode.smart_gap,
            use_log=False,
            quiet=False,
            auxiliary_file=None,
            threads=1
        )
        single_thread_time = time.time() - start_time
        print(f"   Completed in {single_thread_time:.3f} seconds")
    except Exception as e:
        print(f"   Error with single thread: {e}")
        return False
    
    # Test with multiple threads
    print("\n2. Testing with multiple threads (--threads 4):")
    start_time = time.time()
    
    try:
        ck.execute(
            input_file=test_file,
            input_file_format=None,
            output_file="test_output_multi.fa",
            output_file_format=FileFormat.fasta,
            sequence_type=None,
            gaps=0.9,
            gap_characters=None,
            complement=False,
            codon=False,
            ends_only=False,
            mode=TrimmingMode.smart_gap,
            use_log=False,
            quiet=False,
            auxiliary_file=None,
            threads=4
        )
        multi_thread_time = time.time() - start_time
        print(f"   Completed in {multi_thread_time:.3f} seconds")
    except Exception as e:
        print(f"   Error with multiple threads: {e}")
        return False
    
    # Compare outputs to ensure they are identical
    print("\n3. Comparing outputs:")
    with open("test_output_single.fa", "r") as f1, open("test_output_multi.fa", "r") as f2:
        single_content = f1.read()
        multi_content = f2.read()
        
        if single_content == multi_content:
            print("   ✓ Outputs are identical - parallel processing preserves correctness")
        else:
            print("   ✗ Outputs differ - there may be an issue with parallel processing")
            return False
    
    # Clean up test files
    import os
    os.remove("test_output_single.fa")
    os.remove("test_output_multi.fa")
    
    print("\n" + "-" * 50)
    print("Parallel processing test completed successfully!")
    
    if multi_thread_time < single_thread_time:
        speedup = single_thread_time / multi_thread_time
        print(f"Speedup with 4 threads: {speedup:.2f}x")
    
    return True

if __name__ == "__main__":
    success = test_parallel_processing()
    sys.exit(0 if success else 1)
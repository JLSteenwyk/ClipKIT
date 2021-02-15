import pytest
import pytest_mock
from pathlib import Path


import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from clipkit.smart_gap_helper import (
    smart_gap_threshold_determination,
    greatest_diff_in_slopes,
    gap_to_gap_slope, 
    get_gaps_distribution,
    count_and_sort_gaps
)
from clipkit.files import FileFormat

here = Path(__file__)

class TestSmartGapsHelper(object):
    def test_smart_gap_threshold_simple_case(self):
        ## set up
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        expected_gaps = 0.8

        ## execution
        gaps = smart_gap_threshold_determination(alignment)

        ## check results
        assert expected_gaps == gaps

    def test_smart_gap_threshold_standard_case(self):
        ## set up
        alignment = AlignIO.read(f"{here.parent}/examples/EOG091N44M8_aa.fa", "fasta")
        
        ## execution
        gaps = smart_gap_threshold_determination(alignment)
        expected_gaps = 0.8803

        ## check results
        assert expected_gaps == gaps

    def test_greatest_diff_in_slopes_simple_case(self):
        ## set up
        slopes = [0.833333333333333]
        gaps_arr = [[0.8, 1.0], [0.6, 1.0], [0.2, 1.0], [0.0, 3.0]]
        expected_return = 0.8

        ## execution
        greatest_diff = greatest_diff_in_slopes(slopes, gaps_arr)

        ## check results
        assert greatest_diff == expected_return

    def test_greatest_diff_in_slopes_standard_case(self):
        ## set up
        slopes = [
            0.453329706695677,
            0.1903630604288502,
            2.153316106804466,
            0.11466574934067258,
            2.4933133868262236,
            0.6879944960440355,
            0.037924469626292194,
            2.1786492374727793,
            2.4933133868262236,
            4.12796697626421,
            1.6999864001087854,
            1.0319917440660464,
            0.7933269867174482,
            0.11399518940300717,
            0.8026602453847114
        ]
        gaps_arr = [
            [0.9915, 136.0],
            [0.9829, 4.0],
            [0.9573, 5.0],
            [0.9487, 19.0],
            [0.9402, 1.0],
            [0.9316, 22.0],
            [0.9231, 6.0],
            [0.8974, 1.0],
            [0.8889, 19.0],
            [0.8803, 22.0],
            [0.8718, 36.0],
            [0.8632, 15.0],
            [0.8462, 18.0],
            [0.8376, 7.0],
            [0.8205, 2.0],
            [0.812, 7.0],
            [0.8034, 19.0],
            [0.7949, 1.0],
            [0.7692, 288.0],
            [0.7607, 279.0],
            [0.7521, 6.0],
            [0.6667, 5.0],
            [0.1197, 1.0],
            [0.1026, 2.0],
            [0.0855, 1.0],
            [0.0769, 1.0],
            [0.0598, 2.0],
            [0.0342, 1.0],
            [0.0256, 3.0],
            [0.0171, 1.0],
            [0.0085, 1.0],
            [0.0, 95.0]
        ]
        expected_return = 0.8803

        ## execution
        greatest_diff = greatest_diff_in_slopes(slopes, gaps_arr)

        ## check results
        assert greatest_diff == expected_return

    def test_gap_to_gap_slope_simple_case(self):
        ## set up
        gaps_arr = [[0.8, 1.0], [0.6, 1.0], [0.2, 1.0], [0.0, 3.0]]
        alignment_length = 6

        ## execution
        slopes = gap_to_gap_slope(gaps_arr, alignment_length)
        expected_slopes = [0.833333333333333]

        ## check results
        assert expected_slopes == slopes

    def test_get_gaps_distribution(self):
        ## set up
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        alignment_length = 6

        ## execution
        gaps_arr = get_gaps_distribution(alignment, alignment_length)
        expected_gaps_arr = [0.0, 0.6, 0.0, 0.8, 0.0, 0.2]

        ## check results
        assert expected_gaps_arr == gaps_arr

    def test_count_and_sort_gaps(self):
        ## set up
        gaps_dist = [0.0, 0.6, 0.0, 0.8, 0.0, 0.2]

        ## execution
        gaps_arr = count_and_sort_gaps(gaps_dist)
        expected_gaps_arr = [[0.8, 1.0], [0.6, 1.0], [0.2, 1.0], [0.0, 3.0]]

        ## check results
        assert expected_gaps_arr == gaps_arr


        



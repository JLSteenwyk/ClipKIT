import pytest
from pathlib import Path

from Bio import AlignIO
import numpy as np

here = Path(__file__)

from clipkit.modes import TrimmingMode, trim, shouldKeep

class TestModes(object):

    def test_shouldKeep_kpi_gappy_keep(self):
        ## setup
        mode = TrimmingMode.kpi_gappy
        gappyness = 0.00
        gaps = 0.9
        parsimony_informative = True

        assert True == shouldKeep(mode, parsimony_informative, gappyness, gaps)

    def test_shouldKeep_kpi_gappy_trim(self):
        ## setup
        mode = TrimmingMode.kpi_gappy
        gappyness = 0.00
        gaps = 0.9
        parsimony_informative = False

        assert False == shouldKeep(mode, parsimony_informative, gappyness, gaps)

    def test_shouldKeep_gappy_keep(self):
        ## setup
        mode = TrimmingMode.gappy
        gappyness = 0.00
        gaps = 0.9
        parsimony_informative = True

        assert True == shouldKeep(mode, parsimony_informative, gappyness, gaps)

    def test_shouldKeep_gappy_trim(self):
        ## setup
        mode = TrimmingMode.gappy
        gappyness = 0.95
        gaps = 0.9
        parsimony_informative = True

        assert False == shouldKeep(mode, parsimony_informative, gappyness, gaps)

    def test_shouldKeep_kpi_keep(self):
        ## setup
        mode = TrimmingMode.kpi
        gappyness = 0.00
        gaps = 0.9
        parsimony_informative = True

        assert True == shouldKeep(mode, parsimony_informative, gappyness, gaps)

    def test_shouldKeep_kpi_trim(self):
        ## setup
        mode = TrimmingMode.kpi
        gappyness = 0.95
        gaps = 0.9
        parsimony_informative = False

        assert False == shouldKeep(mode, parsimony_informative, gappyness, gaps)

    def test_gappy_mode(self):
        ## setup
        gappyness = 0.00
        parsimony_informative = True
        keepD = {}
        trimD = {}
        i = 2
        gaps = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        use_log = False

        for entry in alignment:
            keepD[entry.id] = np.empty([6], dtype=str)
        trimD = {}
        for entry in alignment:
            trimD[entry.id] = np.empty([6], dtype=str)

        ## execution
        keepD, trimD = trim(gappyness, parsimony_informative, 
            keepD, trimD, i, gaps, alignment, TrimmingMode.gappy, use_log)

        ## check results
        expected_keepD = {
            '1': np.array(['', '', 'G', '', '', '']), 
            '2': np.array(['', '', 'G', '', '', '']),
            '3': np.array(['', '', 'G', '', '', '']),
            '4': np.array(['', '', 'A', '', '', '']),
            '5': np.array(['', '', 'a', '', '', ''])
        }
        expected_trimD = {
            '1': np.array(['', '', '', '', '', '']),
            '2': np.array(['', '', '', '', '', '']),
            '3': np.array(['', '', '', '', '', '']),
            '4': np.array(['', '', '', '', '', '']),
            '5': np.array(['', '', '', '', '', ''])
        }

        assert expected_keepD.keys() == keepD.keys()
        assert all(np.array_equal(expected_keepD[key], keepD[key]) for key in expected_keepD)
        assert expected_trimD.keys() == trimD.keys()
        assert all(np.array_equal(expected_trimD[key], trimD[key]) for key in expected_trimD)

    def test_kpi_gappy_mode(self):
        ## setup
        gappyness = 0.6
        parsimony_informative = False
        keepD = {}
        trimD = {}
        i = 1
        gaps = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        use_log = False

        for entry in alignment:
            keepD[entry.id] = np.empty([6], dtype=str)
        trimD = {}
        for entry in alignment:
            trimD[entry.id] = np.empty([6], dtype=str)

        ## execution
        keepD, trimD = trim(gappyness, parsimony_informative, 
            keepD, trimD, i, gaps, alignment, TrimmingMode.kpi_gappy, use_log)

        ## check results
        expected_keepD = {
            '1': np.array(['', '', '', '', '', '']),
            '2': np.array(['', '', '', '', '', '']),
            '3': np.array(['', '', '', '', '', '']),
            '4': np.array(['', '', '', '', '', '']),
            '5': np.array(['', '', '', '', '', ''])
        }
        expected_trimD = {
            '1': np.array(['', '-', '', '', '', '']),
            '2': np.array(['', '-', '', '', '', '']),
            '3': np.array(['', '-', '', '', '', '']),
            '4': np.array(['', 'G', '', '', '', '']),
            '5': np.array(['', 'C', '', '', '', ''])
        }

        assert expected_keepD.keys() == keepD.keys()
        assert all(np.array_equal(expected_keepD[key], keepD[key]) for key in expected_keepD)
        assert expected_trimD.keys() == trimD.keys()
        assert all(np.array_equal(expected_trimD[key], trimD[key]) for key in expected_trimD)


    def test_kpi_mode(self):
        ## setup
        gappyness = 0.2
        parsimony_informative = True
        keepD = {}
        trimD = {}
        i = 5
        gaps = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        use_log = False

        for entry in alignment:
            keepD[entry.id] = np.empty([6], dtype=str)
        trimD = {}
        for entry in alignment:
            trimD[entry.id] = np.empty([6], dtype=str)

        ## execution
        keepD, trimD = trim(gappyness, parsimony_informative, 
            keepD, trimD, i, gaps, alignment, TrimmingMode.kpi, use_log)

        ## check results
        expected_keepD = {
            '1': np.array(['', '', '', '', '', 'T']),
            '2': np.array(['', '', '', '', '', 'T']),
            '3': np.array(['', '', '', '', '', 'A']),
            '4': np.array(['', '', '', '', '', 'A']),
            '5': np.array(['', '', '', '', '', '-'])
        }
        expected_trimD = {
            '1': np.array(['', '', '', '', '', '']),
            '2': np.array(['', '', '', '', '', '']),
            '3': np.array(['', '', '', '', '', '']),
            '4': np.array(['', '', '', '', '', '']),
            '5': np.array(['', '', '', '', '', ''])
        }

        assert expected_keepD.keys() == keepD.keys()
        assert all(np.array_equal(expected_keepD[key], keepD[key]) for key in expected_keepD)
        assert expected_trimD.keys() == trimD.keys()
        assert all(np.array_equal(expected_trimD[key], trimD[key]) for key in expected_trimD)

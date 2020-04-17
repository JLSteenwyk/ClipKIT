import pytest
from Bio import AlignIO
from pathlib import Path

here = Path(__file__)

from clipkit.modes import gappy_mode
from clipkit.modes import kpi_gappy_mode
from clipkit.modes import kpi_mode

class TestModes(object):

    def test_gappy_mode(self):
        ## setup
        gappyness = 0.00
        parsimony_informative = True
        keepD = {}
        trimD = {}
        i = 2
        gaps = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")

        for entry in alignment:
            keepD[entry.id] = []
        trimD = {}
        for entry in alignment:
            trimD[entry.id] = []

        ## execution
        keepD, trimD = gappy_mode(gappyness, parsimony_informative, 
            keepD, trimD, i, gaps, alignment)

        ## check results
        expected_keepD_keep = {'1': ['G'], '2': ['G'], '3': ['G'], '4': ['A'], '5': ['a']}
        expected_trimD_keep = {'1': [], '2': [], '3': [], '4': [], '5': []}
        expected_logArr = [
            ['3', 'keep', 'PI', '0.0']
            ]
        assert keepD == expected_keepD_keep
        assert trimD == expected_trimD_keep
        # TODO: create logging fixture so we can check expected log output
        # assert logArr == expected_logArr 

    def test_kpi_gappy_mode(self):
        ## setup
        gappyness = 0.6
        parsimony_informative = False
        keepD = {}
        trimD = {}
        i = 1
        gaps = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")

        for entry in alignment:
            keepD[entry.id] = []
        trimD = {}
        for entry in alignment:
            trimD[entry.id] = []

        ## execution
        keepD, trimD = kpi_gappy_mode(gappyness, parsimony_informative, 
            keepD, trimD, i, gaps, alignment)

        ## check results
        expected_keepD = {'1': [], '2': [], '3': [], '4': [], '5': []}
        expected_trimD = {'1': ['-'], '2': ['-'], '3': ['-'], '4': ['G'], '5': ['C']}
        expected_logArr = [
            ['2', 'trim', 'nPI', '0.6']
            ]
        assert keepD == expected_keepD
        assert trimD == expected_trimD
        # TODO: create logging fixture so we can check expected log output
        # assert logArr == expected_logArr 
    
    def test_kpi_mode(self):
        ## setup
        gappyness = 0.2
        parsimony_informative = True
        keepD = {}
        trimD = {}
        i = 5
        gaps = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")

        for entry in alignment:
            keepD[entry.id] = []
        trimD = {}
        for entry in alignment:
            trimD[entry.id] = []

        ## execution
        keepD, trimD = kpi_mode(gappyness, parsimony_informative, 
            keepD, trimD, i, gaps, alignment)

        ## check results
        expected_keepD = {'1': ['T'], '2': ['T'], '3': ['A'], '4': ['A'], '5': ['-']}
        expected_trimD = {'1': [], '2': [], '3': [], '4': [], '5': []}
        expected_logArr = [
            ['6', 'keep', 'PI', '0.2']
            ]
        assert keepD == expected_keepD
        assert trimD == expected_trimD
        # TODO: create logging fixture so we can check expected log output
        # assert logArr == expected_logArr 


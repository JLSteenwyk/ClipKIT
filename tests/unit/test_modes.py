import pytest
from pathlib import Path

import numpy as np
from Bio import AlignIO
from clipkit.modes import SiteClassificationType, TrimmingMode, trim, should_keep_site
from clipkit.msa import MSA
from clipkit.helpers import create_keep_and_trim_msas


here = Path(__file__)


class TestModes(object):
    def test_should_keep_site_kpi_gappy_keep(self):
        ## setup
        mode = TrimmingMode.kpi_gappy
        gappyness = 0.00
        gaps = 0.9

        assert (
            should_keep_site(
                mode, SiteClassificationType.parsimony_informative, gappyness, gaps
            )
            is True
        )

    def test_should_keep_site_kpi_gappy_trim(self):
        ## setup
        mode = TrimmingMode.kpi_gappy
        gappyness = 0.00
        gaps = 0.9

        assert (
            should_keep_site(mode, SiteClassificationType.constant, gappyness, gaps)
            is False
        )

    def test_should_keep_site_gappy_keep(self):
        ## setup
        mode = TrimmingMode.gappy
        site_classification_type = SiteClassificationType.parsimony_informative
        gappyness = 0.00
        gaps_threshold = 0.9

        assert (
            should_keep_site(mode, site_classification_type, gappyness, gaps_threshold)
            is True
        )

    def test_should_keep_site_gappy_trim(self):
        ## setup
        mode = TrimmingMode.gappy
        site_classification_type = SiteClassificationType.parsimony_informative
        gappyness = 0.95
        gaps_threshold = 0.9

        assert (
            should_keep_site(mode, site_classification_type, gappyness, gaps_threshold)
            is False
        )

    def test_should_keep_site_kpi_keep(self):
        ## setup
        mode = TrimmingMode.kpi
        site_classification_type = SiteClassificationType.parsimony_informative
        gappyness = 0.00
        gaps_threshold = 0.9

        assert should_keep_site(
            mode, site_classification_type, gappyness, gaps_threshold
        )

    def test_should_keep_site_kpi_trim(self):
        ## setup
        mode = TrimmingMode.kpi
        site_classification_type = SiteClassificationType.other
        gappyness = 0.95
        gaps_threshold = 0.9

        assert (
            should_keep_site(mode, site_classification_type, gappyness, gaps_threshold)
            is False
        )

    def test_should_keep_site_kpic_keep(self):
        ## setup
        mode = TrimmingMode.kpic
        site_classification_type = SiteClassificationType.constant
        gappyness = 0.95
        gaps_threshold = 0.9

        assert should_keep_site(
            mode, site_classification_type, gappyness, gaps_threshold
        )

    def test_should_keep_site_kpic_trim(self):
        ## setup
        mode = TrimmingMode.kpic
        site_classification_type = SiteClassificationType.other
        gappyness = 0.95
        gaps_threshold = 0.9

        assert (
            should_keep_site(mode, site_classification_type, gappyness, gaps_threshold)
            is False
        )

    def test_should_keep_site_kpic_gappy_keep(self):
        ## setup
        mode = TrimmingMode.kpic_gappy
        site_classification_type = SiteClassificationType.constant
        gappyness = 0.2
        gaps_threshold = 0.9

        assert should_keep_site(
            mode, site_classification_type, gappyness, gaps_threshold
        )

    def test_should_keep_site_kpic_gappy_trim(self):
        ## setup
        mode = TrimmingMode.kpic_gappy
        site_classification_type = SiteClassificationType.constant
        gappyness = 0.95
        gaps_threshold = 0.9

        assert (
            should_keep_site(mode, site_classification_type, gappyness, gaps_threshold)
            is False
        )

    def test_gappy_mode(self):
        ## setup
        gappyness = 0
        site_classification_type = SiteClassificationType.constant
        site_classification_counts = {
            SiteClassificationType.parsimony_informative: 0,
            SiteClassificationType.constant: 0,
            SiteClassificationType.singleton: 0,
            SiteClassificationType.other: 0,
        }
        i = 2
        gaps_threshold = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        use_log = False

        keep_msa, trim_msa = create_keep_and_trim_msas(alignment, True)

        ## execution
        keep_msa, trim_msa = trim(
            gappyness,
            site_classification_type,
            site_classification_counts,
            keep_msa,
            trim_msa,
            i,
            gaps_threshold,
            alignment,
            TrimmingMode.gappy,
            use_log,
        )

        ## check results
        expected_keep_data = {
            "1": np.array([b"", b"", b"G", b"", b"", b""]),
            "2": np.array([b"", b"", b"G", b"", b"", b""]),
            "3": np.array([b"", b"", b"G", b"", b"", b""]),
            "4": np.array([b"", b"", b"A", b"", b"", b""]),
            "5": np.array([b"", b"", b"a", b"", b"", b""]),
        }
        expected_trim_data = {
            "1": np.array([b"", b"", b"", b"", b"", b""]),
            "2": np.array([b"", b"", b"", b"", b"", b""]),
            "3": np.array([b"", b"", b"", b"", b"", b""]),
            "4": np.array([b"", b"", b"", b"", b"", b""]),
            "5": np.array([b"", b"", b"", b"", b"", b""]),
        }

        assert expected_keep_data.keys() == keep_msa._data.keys()
        assert all(
            np.array_equal(expected_keep_data[key], keep_msa._data[key])
            for key in expected_keep_data
        )
        assert expected_trim_data.keys() == trim_msa._data.keys()
        assert all(
            np.array_equal(expected_trim_data[key], trim_msa._data[key])
            for key in expected_trim_data
        )

    def test_kpi_gappy_mode(self):
        ## setup
        gappyness = 0
        site_classification_type = SiteClassificationType.constant
        site_classification_counts = {
            SiteClassificationType.parsimony_informative: 0,
            SiteClassificationType.constant: 0,
            SiteClassificationType.singleton: 0,
            SiteClassificationType.other: 0,
        }
        i = 1
        gaps_threshold = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        use_log = False

        keep_msa, trim_msa = create_keep_and_trim_msas(alignment, True)

        ## execution
        keep_msa, trim_msa = trim(
            gappyness,
            site_classification_type,
            site_classification_counts,
            keep_msa,
            trim_msa,
            i,
            gaps_threshold,
            alignment,
            TrimmingMode.kpi_gappy,
            use_log,
        )

        ## check results
        expected_keep_data = {
            "1": np.array([b"", b"", b"", b"", b"", b""]),
            "2": np.array([b"", b"", b"", b"", b"", b""]),
            "3": np.array([b"", b"", b"", b"", b"", b""]),
            "4": np.array([b"", b"", b"", b"", b"", b""]),
            "5": np.array([b"", b"", b"", b"", b"", b""]),
        }
        expected_trim_data = {
            "1": np.array([b"", b"-", b"", b"", b"", b""]),
            "2": np.array([b"", b"-", b"", b"", b"", b""]),
            "3": np.array([b"", b"-", b"", b"", b"", b""]),
            "4": np.array([b"", b"G", b"", b"", b"", b""]),
            "5": np.array([b"", b"C", b"", b"", b"", b""]),
        }

        assert expected_keep_data.keys() == keep_msa._data.keys()
        assert all(
            np.array_equal(expected_keep_data[key], keep_msa._data[key])
            for key in expected_keep_data
        )
        assert expected_trim_data.keys() == trim_msa._data.keys()
        assert all(
            np.array_equal(expected_trim_data[key], trim_msa._data[key])
            for key in expected_trim_data
        )

    def test_kpi_mode(self):
        ## setup
        gappyness = 0
        site_classification_type = SiteClassificationType.parsimony_informative
        site_classification_counts = {
            SiteClassificationType.parsimony_informative: 0,
            SiteClassificationType.constant: 0,
            SiteClassificationType.singleton: 0,
            SiteClassificationType.other: 0,
        }
        i = 5
        gaps_threshold = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        use_log = False

        keep_msa, trim_msa = create_keep_and_trim_msas(alignment, True)

        ## execution
        keep_msa, trim_msa = trim(
            gappyness,
            site_classification_type,
            site_classification_counts,
            keep_msa,
            trim_msa,
            i,
            gaps_threshold,
            alignment,
            TrimmingMode.kpi,
            use_log,
        )

        ## check results
        expected_keep_data = {
            "1": np.array([b"", b"", b"", b"", b"", b"T"]),
            "2": np.array([b"", b"", b"", b"", b"", b"T"]),
            "3": np.array([b"", b"", b"", b"", b"", b"A"]),
            "4": np.array([b"", b"", b"", b"", b"", b"A"]),
            "5": np.array([b"", b"", b"", b"", b"", b"-"]),
        }
        expected_trim_data = {
            "1": np.array([b"", b"", b"", b"", b"", b""]),
            "2": np.array([b"", b"", b"", b"", b"", b""]),
            "3": np.array([b"", b"", b"", b"", b"", b""]),
            "4": np.array([b"", b"", b"", b"", b"", b""]),
            "5": np.array([b"", b"", b"", b"", b"", b""]),
        }

        assert expected_keep_data.keys() == keep_msa._data.keys()
        assert all(
            np.array_equal(expected_keep_data[key], keep_msa._data[key])
            for key in expected_keep_data
        )
        assert expected_trim_data.keys() == trim_msa._data.keys()
        assert all(
            np.array_equal(expected_trim_data[key], trim_msa._data[key])
            for key in expected_trim_data
        )

    def test_kpic_mode(self):
        ## setup
        gappyness = 0
        site_classification_type = SiteClassificationType.constant
        site_classification_counts = {
            SiteClassificationType.parsimony_informative: 0,
            SiteClassificationType.constant: 0,
            SiteClassificationType.singleton: 0,
            SiteClassificationType.other: 0,
        }
        i = 0
        gaps_threshold = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        use_log = False

        keep_msa, trim_msa = create_keep_and_trim_msas(alignment, True)

        ## execution
        keep_msa, trim_msa = trim(
            gappyness,
            site_classification_type,
            site_classification_counts,
            keep_msa,
            trim_msa,
            i,
            gaps_threshold,
            alignment,
            TrimmingMode.kpic_gappy,
            use_log,
        )

        ## check results
        expected_keep_data = {
            "1": np.array([b"A", b"", b"", b"", b"", b""]),
            "2": np.array([b"A", b"", b"", b"", b"", b""]),
            "3": np.array([b"A", b"", b"", b"", b"", b""]),
            "4": np.array([b"A", b"", b"", b"", b"", b""]),
            "5": np.array([b"A", b"", b"", b"", b"", b""]),
        }
        expected_trim_data = {
            "1": np.array([b"", b"", b"", b"", b"", b""]),
            "2": np.array([b"", b"", b"", b"", b"", b""]),
            "3": np.array([b"", b"", b"", b"", b"", b""]),
            "4": np.array([b"", b"", b"", b"", b"", b""]),
            "5": np.array([b"", b"", b"", b"", b"", b""]),
        }

        assert expected_keep_data.keys() == keep_msa._data.keys()
        assert all(
            np.array_equal(expected_keep_data[key], keep_msa._data[key])
            for key in expected_keep_data
        )
        assert expected_trim_data.keys() == trim_msa._data.keys()
        assert all(
            np.array_equal(expected_trim_data[key], trim_msa._data[key])
            for key in expected_trim_data
        )

    def test_kpic_gappy_mode(self):
        ## setup
        gappyness = 0.2
        site_classification_type = SiteClassificationType.singleton
        site_classification_counts = {
            SiteClassificationType.parsimony_informative: 0,
            SiteClassificationType.constant: 0,
            SiteClassificationType.singleton: 0,
            SiteClassificationType.other: 0,
        }

        i = 3
        gaps_threshold = 0.9
        alignment = AlignIO.read(f"{here.parent}/examples/simple.fa", "fasta")
        use_log = False

        keep_msa, trim_msa = create_keep_and_trim_msas(alignment, True)

        ## execution
        keep_msa, trim_msa = trim(
            gappyness,
            site_classification_type,
            site_classification_counts,
            keep_msa,
            trim_msa,
            i,
            gaps_threshold,
            alignment,
            TrimmingMode.kpic_gappy,
            use_log,
        )

        ## check results
        expected_keep_data = {
            "1": np.array([b"", b"", b"", b"", b"", b""]),
            "2": np.array([b"", b"", b"", b"", b"", b""]),
            "3": np.array([b"", b"", b"", b"", b"", b""]),
            "4": np.array([b"", b"", b"", b"", b"", b""]),
            "5": np.array([b"", b"", b"", b"", b"", b""]),
        }
        expected_trim_data = {
            "1": np.array([b"", b"", b"", b"T", b"", b""]),
            "2": np.array([b"", b"", b"", b"-", b"", b""]),
            "3": np.array([b"", b"", b"", b"-", b"", b""]),
            "4": np.array([b"", b"", b"", b"-", b"", b""]),
            "5": np.array([b"", b"", b"", b"-", b"", b""]),
        }

        assert expected_keep_data.keys() == keep_msa._data.keys()
        assert all(
            np.array_equal(expected_keep_data[key], keep_msa._data[key])
            for key in expected_keep_data
        )
        assert expected_trim_data.keys() == trim_msa._data.keys()
        assert all(
            np.array_equal(expected_trim_data[key], trim_msa._data[key])
            for key in expected_trim_data
        )

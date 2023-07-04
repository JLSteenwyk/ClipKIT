import pytest
import numpy as np

from clipkit.warnings import (
    warn_if_all_sites_were_trimmed,
    warn_if_entry_contains_only_gaps,
)
from clipkit.helpers import SeqType
from clipkit.msa import MSA


class TestWarnings(object):
    def test_all_sites_trimmed(self, mocker):
        mocked_warning = mocker.patch("clipkit.warnings.logger.warning")

        entries = ["some_id"]
        length = 10
        keep_msa = MSA(entries, length)
        print(vars(keep_msa))
        print(keep_msa.entries[0])
        print(keep_msa._data[entries[0]])
        print(keep_msa._data[entries[0]].size)
        warn_if_all_sites_were_trimmed(keep_msa)

        mocked_warning.assert_called_once_with(
            "WARNING: All sites trimmed from alignment. Please use different parameters."
        )

    def test_gaps_only(self, mocker):
        mocked_warning = mocker.patch("clipkit.warnings.logger.warning")

        entries = ["some_id"]
        length = 1
        keep_msa = MSA(entries, length)
        keep_msa.set_entry_sequence_at_position("some_id", 0, "-")
        warn_if_entry_contains_only_gaps(keep_msa, SeqType.aa)

        mocked_warning.assert_called_once_with(
            """WARNING: header id 'some_id' contains only gaps"""
        )

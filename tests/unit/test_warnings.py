import pytest
import numpy as np

from clipkit.warnings import (
    warn_if_all_sites_were_trimmed,
    warn_if_entry_contains_only_gaps,
)
from clipkit.helpers import SeqType
from clipkit.msa import MSA


class TestWarnings(object):
    @pytest.mark.parametrize(
        "header_info, seq_records, should_warn",
        [
            (
                [
                    {"id": "1", "name": "1", "description": "1"},
                    {"id": "2", "name": "2", "description": "2"}
                ],
                np.array([
                    ["", "", "", "", "", ""],
                    ["", "", "", "", "", ""]
                ]),
                True
            ),
            (
                [
                    {"id": "1", "name": "1", "description": "1"},
                    {"id": "2", "name": "2", "description": "2"}
                ],
                np.array([
                    ["A", "-", "G", "T", "A", "T"],
                    ["A", "-", "G", "-", "A", "T"]
                ]),
                False
            ),
        ]
    )
    def test_warn_all_sites_trimmed(self, mocker, header_info, seq_records, should_warn):
        mocked_warning = mocker.patch("clipkit.warnings.logger.warning")

        msa = MSA(header_info, seq_records)
        warn_if_all_sites_were_trimmed(msa)

        if should_warn:
            mocked_warning.assert_called_once_with(
                "WARNING: All sites trimmed from alignment. Please use different parameters."
            )
        else:
            mocked_warning.assert_not_called()


    @pytest.mark.parametrize(
        "header_info, seq_records, gap_only_header_id",
        [
            (
                [
                    {"id": "1", "name": "1", "description": "1"},
                    {"id": "2", "name": "2", "description": "2"}
                ],
                np.array([
                    ["-", "-", "-", "-", "-", "-"],
                    ["A", "G", "G", "T", "A", "C"]
                ]),
                "1"
            ),
            (
                [
                    {"id": "1", "name": "1", "description": "1"},
                    {"id": "2", "name": "2", "description": "2"}
                ],
                np.array([
                    ["A", "G", "G", "T", "A", "C"],
                    ["-", "-", "-", "-", "-", "-"]
                ]),
                "2"
            ),
            (
                [
                    {"id": "1", "name": "1", "description": "1"},
                    {"id": "2", "name": "2", "description": "2"}
                ],
                np.array([
                    ["A", "-", "G", "T", "A", "T"],
                    ["A", "-", "G", "-", "A", "T"]
                ]),
                None
            ),
        ]
    )
    def test_warn_if_entry_contains_only_gaps(self, mocker, header_info, seq_records, gap_only_header_id):
        mocked_warning = mocker.patch("clipkit.warnings.logger.warning")

        msa = MSA(header_info, seq_records)
        warn_if_entry_contains_only_gaps(msa)

        if gap_only_header_id:
            mocked_warning.assert_called_once_with(
                f"""WARNING: header id '{gap_only_header_id}' contains only gaps"""
            )
        else:
            mocked_warning.assert_not_called()

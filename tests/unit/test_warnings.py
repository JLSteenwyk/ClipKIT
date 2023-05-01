import pytest

from clipkit.warnings import warn_if_all_sites_were_trimmed, warn_if_entry_contains_only_gaps


class TestWarnings(object):
    def test_all_sites_trimmed(self, mocker):
        mocked_warning = mocker.patch('clipkit.warnings.logger.warning')
        keepD = dict(some_id="")
        warn_if_all_sites_were_trimmed(keepD)
        mocked_warning.assert_called_once_with('WARNING: All sites trimmed from alignment. Please use different parameters.')

    def test_gaps_only(self, mocker):
        mocked_warning = mocker.patch('clipkit.warnings.logger.warning')
        keepD = dict(some_id="-")
        warn_if_entry_contains_only_gaps(keepD)
        mocked_warning.assert_called_once_with('''WARNING: header id 'some_id' contains only gaps''')

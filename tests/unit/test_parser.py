import pytest

from clipkit.parser import create_parser

@pytest.fixture
def parser():
    return create_parser()


class TestParser(object):
    def test_required_only(self, parser):
        input_path = 'my/input/file.fa'
        parsed = parser.parse_args([input_path])
        assert parsed.input == input_path

    def test_mode(self, parser):
        input_path = 'my/input/file.fa'
        mode = 'gappy'
        parsed = parser.parse_args([input_path, '-m', mode])
        assert parsed.mode == mode
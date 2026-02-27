import pytest

from clipkit.parser import create_parser


@pytest.fixture
def parser():
    return create_parser()


class TestParser(object):
    def test_required_only(self, parser):
        input_path = "my/input/file.fa"
        parsed = parser.parse_args([input_path])
        assert parsed.input == input_path

    def test_mode(self, parser):
        input_path = "my/input/file.fa"
        mode = "gappy"
        parsed = parser.parse_args([input_path, "-m", mode])
        assert parsed.mode == mode

    def test_mode_entropy(self, parser):
        input_path = "my/input/file.fa"
        mode = "entropy"
        parsed = parser.parse_args([input_path, "-m", mode])
        assert parsed.mode == mode

    def test_mode_gappyout(self, parser):
        input_path = "my/input/file.fa"
        mode = "gappyout"
        parsed = parser.parse_args([input_path, "-m", mode])
        assert parsed.mode == mode

    def test_mode_block_gappy(self, parser):
        input_path = "my/input/file.fa"
        mode = "block-gappy"
        parsed = parser.parse_args([input_path, "-m", mode])
        assert parsed.mode == mode

    def test_mode_composition_bias(self, parser):
        input_path = "my/input/file.fa"
        mode = "composition-bias"
        parsed = parser.parse_args([input_path, "-m", mode])
        assert parsed.mode == mode

    def test_mode_heterotachy(self, parser):
        input_path = "my/input/file.fa"
        mode = "heterotachy"
        parsed = parser.parse_args([input_path, "-m", mode])
        assert parsed.mode == mode

    def test_plot_trim_report_with_no_value(self, parser):
        input_path = "my/input/file.fa"
        parsed = parser.parse_args([input_path, "--plot_trim_report"])
        assert parsed.plot_trim_report == ""

    def test_plot_trim_report_with_value(self, parser):
        input_path = "my/input/file.fa"
        parsed = parser.parse_args([input_path, "--plot_trim_report", "report.html"])
        assert parsed.plot_trim_report == "report.html"

import pytest
from argparse import Namespace
from clipkit.modes import TrimmingMode
from clipkit.args_processing import process_args


@pytest.fixture
def args():
    kwargs = dict(
        complementary=False,
        gaps=None,
        input="tests/integration/samples/simple.fa",
        input_file_format=None,
        log=False,
        mode=None,
        output="output/simple",
        output_file_format=None,
    )
    return Namespace(**kwargs)


class TestArgsProcessing(object):
    def test_process_args_input_file_dne(self, args):
        args.input = "some/file/that/doesnt/exist"
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_in_equals_out(self, args):
        args.output = args.input
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_default_mode(self, args):
        res = process_args(args)
        assert res["mode"] == TrimmingMode.smart_gap

    def test_process_args_default_complementary(self, args):
        args.complementary = None
        res = process_args(args)
        assert res["complement"] is False

    def test_process_args_default_gaps(self, args):
        res = process_args(args)
        assert res["gaps"] == 0.9

    def test_process_args_default_use_logs(self, args):
        args.log = None
        res = process_args(args)
        assert res["use_log"] is False

    def test_process_args_default_output_file(self, args):
        args.output = None
        res = process_args(args)
        assert res["output_file"] == f"{args.input}.clipkit"

    def test_process_args_expected_keywords(self, args):
        res = process_args(args)
        expected_keys = [
            "input_file",
            "output_file",
            "input_file_format",
            "output_file_format",
            "complement",
            "gaps",
            "mode",
            "use_log",
        ]
        assert sorted(res.keys()) == sorted(expected_keys)
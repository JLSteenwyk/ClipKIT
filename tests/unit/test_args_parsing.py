from argparse import Namespace
import pytest

from clipkit.args_processing import process_args
from clipkit.helpers import SeqType
from clipkit.modes import TrimmingMode
from clipkit.settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS


@pytest.fixture
def args():
    kwargs = dict(
        complementary=False,
        codon=False,
        gaps=None,
        input="tests/integration/samples/simple.fa",
        input_file_format=None,
        sequence_type=None,
        log=False,
        mode=None,
        output="output/simple",
        output_file_format=None,
        gap_characters=DEFAULT_NT_GAP_CHARS,
        quiet=True,
        auxiliary_file=None,
        ends_only=False,
        dry_run=False,
        validate_only=False,
        report_json=None,
        plot_trim_report=None,
        threads=1,
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

    def test_process_args_default_entropy_threshold(self, args):
        args.mode = TrimmingMode.entropy
        args.gaps = None
        res = process_args(args)
        assert res["gaps"] == 0.8

    def test_process_args_default_use_logs(self, args):
        args.log = None
        res = process_args(args)
        assert res["use_log"] is False

    def test_process_args_default_output_file(self, args):
        args.output = None
        res = process_args(args)
        assert res["output_file"] == f"{args.input}.clipkit"

    def test_process_args_default_sequence_type(self, args):
        args.sequence_type = None
        res = process_args(args)
        assert res["sequence_type"] is None

    def test_process_args_aa_sequence_type(self, args):
        args.sequence_type = "aa"
        res = process_args(args)
        assert res["sequence_type"] == SeqType.aa

    def test_process_args_capital_aa_sequence_type(self, args):
        args.sequence_type = "AA"
        res = process_args(args)
        assert res["sequence_type"] == SeqType.aa

    def test_process_args_gap_characters_nt(self, args):
        res = process_args(args)
        assert res["gap_characters"] == DEFAULT_NT_GAP_CHARS

    def test_process_args_gap_characters_aa(self, args):
        args.gap_characters = DEFAULT_AA_GAP_CHARS
        res = process_args(args)
        assert res["gap_characters"] == DEFAULT_AA_GAP_CHARS

    def test_process_quiet_true(self, args):
        res = process_args(args)
        assert res["quiet"]

    def test_process_quiet_false(self, args):
        args.quiet = False
        res = process_args(args)
        assert res["quiet"] is False

    def test_ends_only_default(self, args):
        res = process_args(args)
        assert res["ends_only"] is False

    def test_ends_only_true(self, args):
        args.ends_only = True
        res = process_args(args)
        assert res["ends_only"] is True

    def test_process_args_expected_keywords(self, args):
        res = process_args(args)
        expected_keys = [
            "input_file",
            "output_file",
            "input_file_format",
            "output_file_format",
            "complement",
            "codon",
            "sequence_type",
            "gaps",
            "mode",
            "use_log",
            "gap_characters",
            "quiet",
            "auxiliary_file",
            "ends_only",
            "dry_run",
            "validate_only",
            "report_json",
            "plot_trim_report",
            "threads",
        ]
        assert sorted(res.keys()) == sorted(expected_keys)

    def test_dry_run_default_false(self, args):
        res = process_args(args)
        assert res["dry_run"] is False

    def test_validate_only_default_false(self, args):
        res = process_args(args)
        assert res["validate_only"] is False

    def test_report_json_default_none(self, args):
        res = process_args(args)
        assert res["report_json"] is None

    def test_report_json_explicit_path(self, args):
        args.report_json = "some/path/report.json"
        res = process_args(args)
        assert res["report_json"] == "some/path/report.json"

    def test_report_json_default_path_when_flag_has_no_value(self, args):
        args.report_json = ""
        res = process_args(args)
        assert res["report_json"] == f"{args.output}.report.json"

    def test_plot_trim_report_default_none(self, args):
        res = process_args(args)
        assert res["plot_trim_report"] is None

    def test_plot_trim_report_explicit_path(self, args):
        args.plot_trim_report = "some/path/trim_report.html"
        res = process_args(args)
        assert res["plot_trim_report"] == "some/path/trim_report.html"

    def test_plot_trim_report_default_path_when_flag_has_no_value(self, args):
        args.plot_trim_report = ""
        res = process_args(args)
        assert res["plot_trim_report"] == f"{args.output}.trim_report.html"

    def test_incompatible_codon_args(self, args):
        args.codon = True
        args.mode = TrimmingMode.c3
        with pytest.raises(SystemExit):
            process_args(args)

    def test_threads_less_than_one_raises(self, args):
        args.threads = 0
        with pytest.raises(SystemExit):
            process_args(args)

    def test_cst_mode_requires_auxiliary_file(self, args):
        args.mode = TrimmingMode.cst
        args.auxiliary_file = None
        with pytest.raises(SystemExit):
            process_args(args)

    def test_cst_mode_requires_existing_auxiliary_file(self, args):
        args.mode = TrimmingMode.cst
        args.auxiliary_file = "does-not-exist.txt"
        with pytest.raises(SystemExit):
            process_args(args)

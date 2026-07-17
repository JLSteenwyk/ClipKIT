import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from clipkit.modes import TrimmingMode

SCRIPT_PATH = Path(__file__).parents[2] / "scripts" / "run_benchmark_smoke.py"
SPEC = importlib.util.spec_from_file_location("run_benchmark_smoke", SCRIPT_PATH)
BENCHMARK = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = BENCHMARK
SPEC.loader.exec_module(BENCHMARK)


def test_comprehensive_suite_covers_every_trimming_mode():
    covered_modes = {
        BENCHMARK.CASES[name].mode
        for name in BENCHMARK.COMPREHENSIVE_CASE_NAMES
        if BENCHMARK.CASES[name].operation == "trim"
    }

    assert covered_modes == {mode.name for mode in TrimmingMode}


def test_comprehensive_suite_covers_entry_points_and_stop_codon_modes():
    cases = [BENCHMARK.CASES[name] for name in BENCHMARK.COMPREHENSIVE_CASE_NAMES]

    assert {case.kind for case in cases} >= {
        "algorithm",
        "execute",
        "cli",
        "api_path",
        "api_raw",
    }
    assert {case.stop_codons for case in cases} >= {"terminal", "internal", "all"}
    assert {case.threads for case in cases} >= {1, 4}


def test_core_suite_covers_every_target_mode_shape_and_entry_point():
    cases = [BENCHMARK.CASES[name] for name in BENCHMARK.CORE_CASE_NAMES]
    algorithm_cases = [case for case in cases if case.kind == "algorithm"]

    for dataset in (
        "large_aa",
        "medium_nt",
        "generated_dense_aa",
        "generated_dense_nt",
    ):
        assert {
            case.mode for case in algorithm_cases if case.dataset == dataset
        } >= set(BENCHMARK.CORE_MODES)

    assert {case.kind for case in cases} >= {
        "algorithm",
        "execute",
        "cli",
        "api_path",
        "api_raw",
    }
    assert {case.threads for case in cases} >= {1, 4}
    assert {
        case.mode
        for case in cases
        if case.kind == "execute" and case.dataset in {"small_aa", "large_aa"}
    } == set(BENCHMARK.CORE_MODES)


def test_sample_summary_rejects_nondeterministic_output():
    case = BENCHMARK.CASES["gappy_small"]
    samples = [
        {
            "runtime_seconds": 1.0,
            "cpu_seconds": 0.9,
            "peak_rss_bytes": 10,
            "output_sha256": "a",
        },
        {
            "runtime_seconds": 1.1,
            "cpu_seconds": 1.0,
            "peak_rss_bytes": 11,
            "output_sha256": "b",
        },
    ]

    with pytest.raises(RuntimeError, match="output changed"):
        BENCHMARK.summarize_samples(case, samples)


def test_sample_summary_rejects_changed_trim_positions():
    case = BENCHMARK.CASES["core_large_aa_kpic_gappy_algorithm"]
    samples = [
        {
            "runtime_seconds": 1.0,
            "cpu_seconds": 0.9,
            "peak_rss_bytes": 10,
            "output_sha256": "same",
            "trim_positions_sha256": "a",
        },
        {
            "runtime_seconds": 1.1,
            "cpu_seconds": 1.0,
            "peak_rss_bytes": 11,
            "output_sha256": "same",
            "trim_positions_sha256": "b",
        },
    ]

    with pytest.raises(RuntimeError, match="trim_positions_sha256 changed"):
        BENCHMARK.summarize_samples(case, samples)


def test_sample_summary_reports_median_range_and_memory():
    case = BENCHMARK.CASES["gappy_small"]
    samples = [
        {
            "runtime_seconds": 3.0,
            "cpu_seconds": 2.5,
            "peak_rss_bytes": 30,
            "output_sha256": "same",
        },
        {
            "runtime_seconds": 1.0,
            "cpu_seconds": 0.5,
            "peak_rss_bytes": 10,
            "output_sha256": "same",
        },
        {
            "runtime_seconds": 2.0,
            "cpu_seconds": 1.5,
            "peak_rss_bytes": 20,
            "output_sha256": "same",
        },
    ]

    result = BENCHMARK.summarize_samples(case, samples)

    assert result["runtime_seconds_median"] == 2.0
    assert result["runtime_seconds_range"] == 2.0
    assert result["cpu_seconds_median"] == 1.5
    assert result["cpu_seconds_range"] == 2.0
    assert result["peak_rss_bytes_median"] == 20


def _case_input(case_name):
    case = BENCHMARK.CASES[case_name]
    return BENCHMARK.ROOT / BENCHMARK.FILES[case.dataset]


@pytest.mark.parametrize(
    "case_name",
    [
        "construct_small_aa",
        "gappy_small",
        "api_path_gappy_small",
        "api_raw_gappy_small",
        "cli_gappy_small",
        "output_ecomp",
        "input_ecomp",
    ],
)
def test_worker_entry_points_produce_deterministic_metadata(case_name, capsys):
    result = BENCHMARK.run_worker(case_name, _case_input(case_name))
    payload = json.loads(capsys.readouterr().out)

    assert result == 0
    assert payload["runtime_seconds"] >= 0
    assert payload["cpu_seconds"] is None or payload["cpu_seconds"] >= 0
    assert len(payload["output_sha256"]) == 64
    assert payload["effective_threads"] is None or payload["effective_threads"] >= 1


def test_algorithm_worker_covers_smart_gap_classification_metadata(tmp_path, capsys):
    input_path = tmp_path / "alignment.fa"
    input_path.write_text(">a\nAA-C\n>b\nAC-C\n>c\nCAAC\n>d\nCCAC\n")
    case = BENCHMARK.BenchmarkCase(
        "test_kpic_smart_gap",
        "algorithm",
        "temporary",
        mode="kpic_smart_gap",
    )
    BENCHMARK.CASES[case.name] = case
    try:
        assert BENCHMARK.run_worker(case.name, input_path) == 0
    finally:
        BENCHMARK.CASES.pop(case.name)

    payload = json.loads(capsys.readouterr().out)
    assert payload["threshold"] is not None
    assert len(payload["classification_sha256"]) == 64
    assert len(payload["gappyness_sha256"]) == 64


def test_algorithm_worker_masks_requested_stop_codons(tmp_path, capsys):
    input_path = tmp_path / "stops.fa"
    input_path.write_text(">stop\nATGTAATGA\n>control\nATGCAACAA\n")
    case = BENCHMARK.BenchmarkCase(
        "test_stop_codons",
        "algorithm",
        "temporary",
        operation="stop_codons",
        sequence_type="nt",
        codon=True,
        stop_codons="all",
    )
    BENCHMARK.CASES[case.name] = case
    try:
        assert BENCHMARK.run_worker(case.name, input_path) == 0
    finally:
        BENCHMARK.CASES.pop(case.name)

    payload = json.loads(capsys.readouterr().out)
    assert payload["threshold"] is None
    assert len(payload["output_sha256"]) == 64


def test_cli_worker_forwards_codon_and_stop_codon_options(tmp_path):
    input_path = tmp_path / "stops.fa"
    input_path.write_text(">stop\nATGTAATGA\n>control\nATGCAACAA\n")
    case = BENCHMARK.BenchmarkCase(
        "cli_stop_codons",
        "cli",
        "temporary",
        mode="gappy",
        gaps=0.9,
        sequence_type="nt",
        codon=True,
        stop_codons="all",
        ends_only=True,
        input_format="fasta",
        output_format="fasta",
    )

    result = BENCHMARK._cli_worker(case, input_path)

    assert result["comparison_type"] == "exact_bytes"
    assert result["output_bytes"] > 0
    assert len(result["output_sha256"]) == 64


def test_cli_worker_forwards_auxiliary_file():
    case = BENCHMARK.CASES["cst_tiny"]

    result = BENCHMARK._cli_worker(case, _case_input("cst_tiny"))

    assert result["comparison_type"] == "exact_bytes"
    assert result["output_bytes"] > 0


def test_cli_worker_reports_subprocess_failure(monkeypatch):
    monkeypatch.setattr(
        BENCHMARK.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=2, stderr="bad option"),
    )

    with pytest.raises(RuntimeError, match="(?s)ClipKIT CLI failed.*bad option"):
        BENCHMARK._cli_worker(
            BENCHMARK.CASES["cli_gappy_small"],
            _case_input("cli_gappy_small"),
        )


def test_fasta_dimensions_support_multiline_records_and_other_formats(tmp_path):
    input_path = tmp_path / "multiline.fa"
    input_path.write_text(">first\nAC\nGT\n>second\nA-\nGT\n")

    assert BENCHMARK._fasta_dimensions(input_path, None) == (2, 4)
    assert BENCHMARK._fasta_dimensions(input_path, "ecomp") is None


def test_resource_measurements_are_non_negative():
    peak_rss = BENCHMARK._peak_rss_bytes()
    child_peak_rss = BENCHMARK._peak_rss_bytes(children=True)
    child_cpu = BENCHMARK._child_cpu_seconds()

    assert peak_rss is None or peak_rss >= 0
    assert child_peak_rss is None or child_peak_rss >= 0
    assert child_cpu is None or child_cpu >= 0


def test_unknown_worker_kind_is_rejected(tmp_path):
    input_path = tmp_path / "alignment.fa"
    input_path.write_text(">a\nAAAA\n>b\nAAAA\n")
    case = BENCHMARK.BenchmarkCase("unknown_worker_kind", "unsupported", "temporary")
    BENCHMARK.CASES[case.name] = case
    try:
        with pytest.raises(ValueError, match="Unknown benchmark kind"):
            BENCHMARK.run_worker(case.name, input_path)
    finally:
        BENCHMARK.CASES.pop(case.name)

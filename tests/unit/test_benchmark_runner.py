import importlib.util
import sys
from pathlib import Path

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


def test_sample_summary_rejects_nondeterministic_output():
    case = BENCHMARK.CASES["gappy_small"]
    samples = [
        {"runtime_seconds": 1.0, "peak_rss_bytes": 10, "output_sha256": "a"},
        {"runtime_seconds": 1.1, "peak_rss_bytes": 11, "output_sha256": "b"},
    ]

    with pytest.raises(RuntimeError, match="output changed"):
        BENCHMARK.summarize_samples(case, samples)


def test_sample_summary_reports_median_range_and_memory():
    case = BENCHMARK.CASES["gappy_small"]
    samples = [
        {"runtime_seconds": 3.0, "peak_rss_bytes": 30, "output_sha256": "same"},
        {"runtime_seconds": 1.0, "peak_rss_bytes": 10, "output_sha256": "same"},
        {"runtime_seconds": 2.0, "peak_rss_bytes": 20, "output_sha256": "same"},
    ]

    result = BENCHMARK.summarize_samples(case, samples)

    assert result["runtime_seconds_median"] == 2.0
    assert result["runtime_seconds_range"] == 2.0
    assert result["peak_rss_bytes_median"] == 20

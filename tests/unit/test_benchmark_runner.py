import builtins
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


def _benchmark_sample(**overrides):
    sample = {
        "runtime_seconds": 1.0,
        "cpu_seconds": 0.5,
        "peak_rss_bytes": 100,
        "output_sha256": "same-output",
    }
    sample.update(overrides)
    return sample


def test_sample_summary_preserves_optional_metadata_and_missing_metrics():
    case = BENCHMARK.CASES["gappy_small"]
    result = BENCHMARK.summarize_samples(
        case,
        [
            _benchmark_sample(
                cpu_seconds=None,
                peak_rss_bytes=None,
                threshold=0.9,
                effective_threads=1,
            )
        ],
    )

    assert result["cpu_seconds_median"] is None
    assert result["peak_rss_bytes_median"] is None
    assert result["threshold"] == 0.9
    assert result["effective_threads"] == 1


def test_invoke_worker_runs_fresh_process():
    case = BENCHMARK.CASES["gappy_small"]

    result = BENCHMARK._invoke_worker(
        BENCHMARK.ROOT,
        case,
        _case_input(case.name),
    )

    assert result["runtime_seconds"] >= 0
    assert len(result["output_sha256"]) == 64


def test_invoke_worker_reports_process_failure(monkeypatch):
    monkeypatch.setattr(
        BENCHMARK.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=1, stderr="worker failed"),
    )
    case = BENCHMARK.CASES["gappy_small"]

    with pytest.raises(
        RuntimeError, match="(?s)Benchmark worker failed.*worker failed"
    ):
        BENCHMARK._invoke_worker(
            BENCHMARK.ROOT,
            case,
            _case_input(case.name),
        )


def test_run_case_discards_warmups_and_summarizes_repetitions(monkeypatch):
    calls = []

    def invoke_worker(source_root, case, input_path):
        calls.append((source_root, case.name, input_path))
        return _benchmark_sample()

    monkeypatch.setattr(BENCHMARK, "_invoke_worker", invoke_worker)
    case = BENCHMARK.CASES["gappy_small"]

    result = BENCHMARK.run_case(
        BENCHMARK.ROOT,
        case,
        _case_input(case.name),
        warmups=2,
        repetitions=3,
    )

    assert len(calls) == 5
    assert result["repetitions"] == 3
    assert result["runtime_seconds_median"] == 1.0


def test_generated_fixtures_are_deterministic_and_well_formed(tmp_path, monkeypatch):
    def shortened_range(*args):
        if args == (384,) or args == (512,):
            return builtins.range(2)
        if args == (5000,) or args == (9000,):
            return builtins.range(20)
        return builtins.range(*args)

    monkeypatch.setattr(BENCHMARK, "range", shortened_range, raising=False)

    generated = BENCHMARK._write_generated_fixtures(tmp_path)

    assert set(generated) == {
        "generated_dense_aa",
        "generated_dense_nt",
        "generated_sparse_aa",
        "generated_stop_nt",
    }
    assert BENCHMARK._fasta_dimensions(generated["generated_dense_aa"], None) == (
        2,
        20,
    )
    assert BENCHMARK._fasta_dimensions(generated["generated_dense_nt"], None) == (
        2,
        20,
    )
    assert BENCHMARK._fasta_dimensions(generated["generated_sparse_aa"], None) == (
        2,
        20,
    )
    assert BENCHMARK._fasta_dimensions(generated["generated_stop_nt"], None) == (
        2,
        6000,
    )
    assert "---" in generated["generated_stop_nt"].read_text()


def test_input_paths_resolve_static_and_generated_datasets(tmp_path):
    generated = {
        name: tmp_path / f"{name}.fa"
        for name, relative in BENCHMARK.FILES.items()
        if relative is None
    }

    paths = BENCHMARK._input_paths(BENCHMARK.ROOT, generated)

    assert paths["small_aa"] == BENCHMARK.ROOT / BENCHMARK.FILES["small_aa"]
    assert paths["generated_dense_aa"] == generated["generated_dense_aa"]


def test_git_revision_reports_repository_and_non_repository(tmp_path, monkeypatch):
    revision = BENCHMARK._git_revision(BENCHMARK.ROOT)

    assert revision is not None
    assert len(revision) == 40

    monkeypatch.setattr(
        BENCHMARK.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=128, stdout=""),
    )
    assert BENCHMARK._git_revision(tmp_path) is None


def _generated_paths(tmp_path):
    return {
        name: tmp_path / f"{name}.fa"
        for name, relative in BENCHMARK.FILES.items()
        if relative is None
    }


def test_run_source_reports_each_case_and_revision(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        BENCHMARK,
        "run_case",
        lambda *args, **kwargs: {
            "runtime_seconds_median": 0.25,
            "output_sha256": "equivalent",
        },
    )
    monkeypatch.setattr(BENCHMARK, "_git_revision", lambda source_root: "revision")

    result = BENCHMARK._run_source(
        BENCHMARK.ROOT,
        ("gappy_small", "cli_gappy_small"),
        _generated_paths(tmp_path),
        warmups=0,
        repetitions=1,
    )

    assert result["git_revision"] == "revision"
    assert set(result["cases"]) == {"gappy_small", "cli_gappy_small"}
    output = capsys.readouterr().out
    assert "gappy_small: 0.250000s median" in output
    assert "cli_gappy_small: 0.250000s median" in output


def test_equivalence_groups_reject_different_outputs():
    results = {
        "gappy_small": {"output_sha256": "first"},
        "cli_gappy_small": {"output_sha256": "second"},
    }

    with pytest.raises(RuntimeError, match="Equivalent cases produced different"):
        BENCHMARK._verify_equivalence_groups(results)


def _comparison_case(**overrides):
    result = {
        "output_sha256": "same",
        "output_bytes": 100,
        "threshold": 0.9,
        "comparison_type": "exact_bytes",
        "keep_positions_sha256": "keep",
        "trim_positions_sha256": "trim",
        "classification_sha256": "classification",
        "gappyness_sha256": "gappyness",
        "effective_threads": 1,
        "runtime_seconds_median": 2.0,
        "peak_rss_bytes_median": 100,
        "cpu_seconds_median": 1.0,
    }
    result.update(overrides)
    return result


def test_compare_sources_reports_performance_changes():
    candidate = {"cases": {"case": _comparison_case()}}
    reference = {
        "cases": {
            "case": _comparison_case(
                runtime_seconds_median=4.0,
                peak_rss_bytes_median=200,
                cpu_seconds_median=2.0,
            )
        }
    }

    result = BENCHMARK.compare_sources(candidate, reference)["case"]

    assert result == {
        "speedup": 2.0,
        "runtime_change_percent": -50.0,
        "cpu_change_percent": -50.0,
        "peak_rss_change_percent": -50.0,
    }


def test_compare_sources_supports_unavailable_resource_metrics():
    candidate = {
        "cases": {
            "case": _comparison_case(
                cpu_seconds_median=None,
                peak_rss_bytes_median=None,
            )
        }
    }
    reference = {
        "cases": {
            "case": _comparison_case(
                cpu_seconds_median=None,
                peak_rss_bytes_median=None,
            )
        }
    }

    result = BENCHMARK.compare_sources(candidate, reference)["case"]

    assert result["cpu_change_percent"] is None
    assert result["peak_rss_change_percent"] is None


def test_compare_sources_rejects_biological_differences():
    candidate = {"cases": {"case": _comparison_case(output_sha256="changed")}}
    reference = {"cases": {"case": _comparison_case()}}

    with pytest.raises(RuntimeError, match="differs.*output_sha256"):
        BENCHMARK.compare_sources(candidate, reference)


def test_dependency_versions_report_missing_packages(monkeypatch):
    def version(name):
        if name == "numpy":
            return "1.2.3"
        raise BENCHMARK.importlib.metadata.PackageNotFoundError(name)

    monkeypatch.setattr(BENCHMARK.importlib.metadata, "version", version)

    assert BENCHMARK._dependency_versions() == {
        "numpy": "1.2.3",
        "biopython": None,
    }


@pytest.mark.parametrize(
    "arguments, message",
    [
        (["--worker"], "--worker requires --case and --input"),
        (["--repetitions", "0"], "--repetitions must be at least 1"),
        (["--warmups", "-1"], "--warmups cannot be negative"),
    ],
)
def test_main_rejects_invalid_worker_and_sampling_arguments(
    arguments, message, monkeypatch, capsys
):
    monkeypatch.setattr(sys, "argv", [str(SCRIPT_PATH), *arguments])

    with pytest.raises(SystemExit) as error:
        BENCHMARK.main()

    assert error.value.code == 2
    assert message in capsys.readouterr().err


def test_main_dispatches_worker_mode(monkeypatch, tmp_path):
    input_path = tmp_path / "input.fa"
    input_path.write_text(">a\nAAAA\n>b\nAAAA\n")
    received = {}

    def run_worker(case_name, worker_input):
        received["case_name"] = case_name
        received["input"] = worker_input
        return 7

    monkeypatch.setattr(BENCHMARK, "run_worker", run_worker)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT_PATH),
            "--worker",
            "--case",
            "gappy_small",
            "--input",
            str(input_path),
        ],
    )

    assert BENCHMARK.main() == 7
    assert received == {"case_name": "gappy_small", "input": input_path}


@pytest.mark.parametrize(
    "suite, expected_cases",
    [
        ("smoke", BENCHMARK.SMOKE_CASE_NAMES),
        ("full", BENCHMARK.FULL_CASE_NAMES),
        ("core", BENCHMARK.CORE_CASE_NAMES),
        ("comprehensive", BENCHMARK.COMPREHENSIVE_CASE_NAMES),
    ],
)
def test_main_writes_report_for_each_suite(
    suite, expected_cases, tmp_path, monkeypatch
):
    output_path = tmp_path / f"{suite}.json"
    received = {}

    def run_source(source_root, case_names, generated, warmups, repetitions):
        received["source_root"] = source_root
        received["case_names"] = case_names
        received["warmups"] = warmups
        received["repetitions"] = repetitions
        return {"source_root": str(source_root), "git_revision": None, "cases": {}}

    monkeypatch.setattr(BENCHMARK, "_write_generated_fixtures", lambda path: {})
    monkeypatch.setattr(BENCHMARK, "_run_source", run_source)
    monkeypatch.setattr(BENCHMARK, "_dependency_versions", lambda: {"numpy": "test"})
    monkeypatch.setattr(BENCHMARK.time, "time", lambda: 1234567890)
    monkeypatch.setattr(BENCHMARK.platform, "platform", lambda: "test-platform")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT_PATH),
            "--output",
            str(output_path),
            "--suite",
            suite,
            "--warmups",
            "0",
            "--repetitions",
            "1",
        ],
    )

    assert BENCHMARK.main() == 0
    payload = json.loads(output_path.read_text())
    assert received["source_root"] == BENCHMARK.ROOT.resolve()
    assert received["case_names"] == expected_cases
    assert received["warmups"] == 0
    assert received["repetitions"] == 1
    assert payload["timestamp_unix"] == 1234567890
    assert payload["suite"] == suite
    assert payload["platform"] == "test-platform"
    assert payload["dependencies"] == {"numpy": "test"}
    assert payload["reference"] is None
    assert payload["candidate_vs_reference"] is None


def test_main_compares_candidate_and_reference_sources(tmp_path, monkeypatch):
    output_path = tmp_path / "comparison.json"
    candidate_root = tmp_path / "candidate"
    reference_root = tmp_path / "reference"
    candidate_root.mkdir()
    reference_root.mkdir()
    source_results = [
        {
            "source_root": str(candidate_root),
            "git_revision": "candidate",
            "cases": {"case": _comparison_case()},
        },
        {
            "source_root": str(reference_root),
            "git_revision": "reference",
            "cases": {
                "case": _comparison_case(runtime_seconds_median=4.0),
            },
        },
    ]

    monkeypatch.setattr(BENCHMARK, "_write_generated_fixtures", lambda path: {})
    monkeypatch.setattr(
        BENCHMARK,
        "_run_source",
        lambda *args, **kwargs: source_results.pop(0),
    )
    monkeypatch.setattr(BENCHMARK, "_dependency_versions", lambda: {})
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT_PATH),
            "--output",
            str(output_path),
            "--source-root",
            str(candidate_root),
            "--compare-root",
            str(reference_root),
            "--warmups",
            "0",
            "--repetitions",
            "1",
        ],
    )

    assert BENCHMARK.main() == 0
    payload = json.loads(output_path.read_text())
    assert payload["reference"]["git_revision"] == "reference"
    assert payload["candidate_vs_reference"]["case"]["speedup"] == 2.0

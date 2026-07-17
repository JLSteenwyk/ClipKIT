#!/usr/bin/env python3
"""Reproducible ClipKIT correctness and performance benchmark runner.

Each measured sample runs in a fresh Python process.  Warm-up samples are
discarded, and the retained samples must produce byte-for-byte identical
results.  Passing ``--compare-root`` additionally verifies a candidate source
tree against a reference source tree before reporting relative performance.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import os
import platform
import random
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]

FILES = {
    "tiny_nt": "tests/integration/samples/simple.fa",
    "small_aa": "tests/integration/samples/EOG091N44M8_aa.fa",
    "medium_aa": "tests/integration/samples/EOG092C4VOX_aa_aln.fasta",
    "large_aa": "tests/integration/samples/EOG092C0CZK_aa_aln.fasta",
    "medium_nt": "tests/integration/samples/EOG091N44M8_nt.fa",
    "codon_nt": "tests/integration/samples/12_YIL115C_Anc_2.253_codon_aln.fasta",
    "ecomp_aa": "tests/integration/samples/12_YIL115C_Anc_2.253_aa_aln.ecomp",
    "generated_sparse_aa": None,
    "generated_stop_nt": None,
}


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    kind: str
    dataset: str
    operation: str = "trim"
    mode: str = "gappy"
    gaps: float | None = 0.9
    sequence_type: str | None = None
    codon: bool = False
    stop_codons: str | None = None
    ends_only: bool = False
    input_format: str | None = None
    output_format: str | None = "fasta"
    auxiliary_file: str | None = None
    threads: int = 1


CASES = {
    case.name: case
    for case in [
        # Construction and the principal algorithm-level hot paths.
        BenchmarkCase("construct_small_aa", "algorithm", "small_aa", "construct"),
        BenchmarkCase("construct_medium_aa", "algorithm", "medium_aa", "construct"),
        BenchmarkCase("construct_large_aa", "algorithm", "large_aa", "construct"),
        BenchmarkCase("construct_medium_nt", "algorithm", "medium_nt", "construct"),
        BenchmarkCase("gappy_large_algorithm", "algorithm", "large_aa", mode="gappy"),
        BenchmarkCase(
            "smart_gap_large_algorithm", "algorithm", "large_aa", mode="smart_gap"
        ),
        BenchmarkCase("kpic_large_algorithm", "algorithm", "large_aa", mode="kpic"),
        BenchmarkCase(
            "entropy_large_algorithm", "algorithm", "large_aa", mode="entropy"
        ),
        BenchmarkCase(
            "composition_bias_large_algorithm",
            "algorithm",
            "large_aa",
            mode="composition_bias",
        ),
        BenchmarkCase(
            "sparse_gappy_algorithm",
            "algorithm",
            "generated_sparse_aa",
            mode="gappy",
        ),
        BenchmarkCase(
            "stop_terminal_algorithm",
            "algorithm",
            "generated_stop_nt",
            "stop_codons",
            sequence_type="nt",
            codon=True,
            stop_codons="terminal",
        ),
        BenchmarkCase(
            "stop_internal_algorithm",
            "algorithm",
            "generated_stop_nt",
            "stop_codons",
            sequence_type="nt",
            codon=True,
            stop_codons="internal",
        ),
        BenchmarkCase(
            "stop_all_algorithm",
            "algorithm",
            "generated_stop_nt",
            "stop_codons",
            sequence_type="nt",
            codon=True,
            stop_codons="all",
        ),
        # End-to-end representative sizes and every trimming mode.
        BenchmarkCase("gappy_small", "execute", "small_aa", mode="gappy"),
        BenchmarkCase("smart_gap_medium", "execute", "medium_aa", mode="smart_gap"),
        BenchmarkCase("gappy_large", "execute", "large_aa", mode="gappy"),
        BenchmarkCase("kpic_large", "execute", "large_aa", mode="kpic"),
        BenchmarkCase("block_gappy_small", "execute", "small_aa", mode="block_gappy"),
        BenchmarkCase("gappyout_small", "execute", "small_aa", mode="gappyout"),
        BenchmarkCase("entropy_small", "execute", "small_aa", mode="entropy", gaps=0.8),
        BenchmarkCase(
            "composition_bias_small",
            "execute",
            "small_aa",
            mode="composition_bias",
            gaps=0.8,
        ),
        BenchmarkCase(
            "heterotachy_tiny",
            "execute",
            "tiny_nt",
            mode="heterotachy",
            gaps=0.8,
            sequence_type="nt",
        ),
        BenchmarkCase("kpi_small", "execute", "small_aa", mode="kpi"),
        BenchmarkCase("kpi_gappy_small", "execute", "small_aa", mode="kpi_gappy"),
        BenchmarkCase(
            "kpi_smart_gap_small", "execute", "small_aa", mode="kpi_smart_gap"
        ),
        BenchmarkCase("kpic_small", "execute", "small_aa", mode="kpic"),
        BenchmarkCase("kpic_gappy_small", "execute", "small_aa", mode="kpic_gappy"),
        BenchmarkCase(
            "kpic_smart_gap_small", "execute", "small_aa", mode="kpic_smart_gap"
        ),
        BenchmarkCase(
            "cst_tiny",
            "execute",
            "tiny_nt",
            mode="cst",
            sequence_type="nt",
            auxiliary_file="tests/integration/samples/cst_ex0.txt",
        ),
        BenchmarkCase(
            "c3_medium_nt",
            "execute",
            "medium_nt",
            mode="c3",
            sequence_type="nt",
        ),
        BenchmarkCase(
            "codon_medium_nt",
            "execute",
            "codon_nt",
            mode="gappy",
            sequence_type="nt",
            codon=True,
        ),
        BenchmarkCase(
            "ends_only_medium_nt",
            "execute",
            "medium_nt",
            mode="gappy",
            sequence_type="nt",
            ends_only=True,
        ),
        BenchmarkCase(
            "threads_1_medium_aa",
            "execute",
            "medium_aa",
            mode="smart_gap",
            threads=1,
        ),
        BenchmarkCase(
            "threads_4_medium_aa",
            "execute",
            "medium_aa",
            mode="smart_gap",
            threads=4,
        ),
        # Public CLI/API paths and stop-codon behavior.
        BenchmarkCase("cli_gappy_small", "cli", "small_aa", mode="gappy"),
        BenchmarkCase("api_path_gappy_small", "api_path", "small_aa", mode="gappy"),
        BenchmarkCase("api_raw_gappy_small", "api_raw", "small_aa", mode="gappy"),
        BenchmarkCase(
            "stop_terminal_end_to_end",
            "execute",
            "generated_stop_nt",
            mode="gappy",
            sequence_type="nt",
            codon=True,
            stop_codons="terminal",
        ),
        BenchmarkCase(
            "stop_internal_end_to_end",
            "execute",
            "generated_stop_nt",
            mode="gappy",
            sequence_type="nt",
            codon=True,
            stop_codons="internal",
        ),
        BenchmarkCase(
            "stop_all_end_to_end",
            "execute",
            "generated_stop_nt",
            mode="gappy",
            sequence_type="nt",
            codon=True,
            stop_codons="all",
        ),
        # Supported input/output formats use a tiny fixture so these measure
        # compatibility without dominating the suite.
        BenchmarkCase(
            "input_ecomp", "execute", "ecomp_aa", mode="gappy", input_format="ecomp"
        ),
        *[
            BenchmarkCase(
                f"output_{format_name}",
                "execute",
                "tiny_nt",
                mode="gappy",
                sequence_type="nt",
                input_format="fasta",
                output_format=format_name,
            )
            for format_name in (
                "fasta",
                "clustal",
                "maf",
                "mauve",
                "phylip",
                "phylip_sequential",
                "phylip_relaxed",
                "stockholm",
                "ecomp",
            )
        ],
    ]
}

SMOKE_CASE_NAMES = ("gappy_small", "smart_gap_medium", "kpic_large")

FULL_CASE_NAMES = (
    "construct_small_aa",
    "construct_medium_aa",
    "construct_large_aa",
    "gappy_large_algorithm",
    "smart_gap_large_algorithm",
    "kpic_large_algorithm",
    "entropy_large_algorithm",
    "composition_bias_large_algorithm",
    "gappy_small",
    "smart_gap_medium",
    "gappy_large",
    "kpic_large",
    "c3_medium_nt",
)

COMPREHENSIVE_CASE_NAMES = tuple(CASES)

# Cases within each group must have identical biological output.  This catches
# nondeterminism or behavior changes caused by thread selection and entry point.
EQUIVALENCE_GROUPS = (
    ("smart_gap_medium", "threads_1_medium_aa", "threads_4_medium_aa"),
    ("gappy_small", "cli_gappy_small"),
)


def _peak_rss_bytes(children: bool = False) -> int | None:
    try:
        import resource
    except ImportError:
        return None

    target = resource.RUSAGE_CHILDREN if children else resource.RUSAGE_SELF
    peak_rss = resource.getrusage(target).ru_maxrss
    return peak_rss if sys.platform == "darwin" else peak_rss * 1024


def _hash_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _msa_digest(msa, extra: bytes = b"") -> str:
    payload = (
        msa.seq_records.tobytes()
        + msa._site_positions_to_keep.astype("int64").tobytes()
        + msa._site_positions_to_trim.astype("int64").tobytes()
        + extra
    )
    return _hash_bytes(payload)


def _api_digest(trim_run, stats) -> str:
    payload = {
        "headers": trim_run.msa.header_info,
        "sequences": [str(record.seq) for record in trim_run.trimmed],
        "sites_kept": [int(value) for value in trim_run.msa._site_positions_to_keep],
        "sites_trimmed": [int(value) for value in trim_run.msa._site_positions_to_trim],
        "stats": stats.summary,
        "stop_codons": trim_run.stop_codon_masking.summary,
    }
    return _hash_bytes(json.dumps(payload, sort_keys=True).encode("utf-8"))


def _output_digest(output_path: Path, output_format: str | None) -> tuple[str, str]:
    if output_format != "ecomp":
        return _hash_bytes(output_path.read_bytes()), "exact_bytes"

    # ECOMP's gzip fallback carries a creation timestamp, so two archives made
    # from the same alignment need not be byte-identical.  Compare their fully
    # decoded biological content and metadata instead.
    from clipkit.ecomp import read_ecomp

    alignment, metadata = read_ecomp(output_path)
    payload = {
        "records": [
            (record.id, record.description, str(record.seq)) for record in alignment
        ],
        "metadata": metadata,
    }
    canonical = json.dumps(payload, default=str, sort_keys=True).encode("utf-8")
    return _hash_bytes(canonical), "decoded_ecomp"


def _read_alignment(input_path: Path, input_format: str | None):
    if input_format == "ecomp":
        from clipkit.ecomp import read_ecomp

        return read_ecomp(input_path)[0]

    from Bio import AlignIO

    return AlignIO.read(input_path, input_format or "fasta")


def _algorithm_worker(case: BenchmarkCase, input_path: Path) -> dict[str, Any]:
    from clipkit.modes import StopCodonMode, TrimmingMode
    from clipkit.msa import MSA
    from clipkit.settings import DEFAULT_AA_GAP_CHARS, DEFAULT_NT_GAP_CHARS
    from clipkit.smart_gap_helper import smart_gap_threshold_determination

    alignment = _read_alignment(input_path, case.input_format)
    gaps = DEFAULT_NT_GAP_CHARS if case.sequence_type == "nt" else DEFAULT_AA_GAP_CHARS
    start = time.perf_counter()
    msa = MSA.from_bio_msa(alignment, gaps)
    threshold = None
    extra = b""

    if case.operation == "stop_codons":
        stats = msa.mask_stop_codons(StopCodonMode(case.stop_codons))
        extra = json.dumps(stats.summary, sort_keys=True).encode("utf-8")
    elif case.operation == "trim":
        mode = TrimmingMode[case.mode]
        threshold = case.gaps
        if case.mode in {"smart_gap", "kpi_smart_gap", "kpic_smart_gap"}:
            kwargs = {"seq_records": msa.seq_records}
            if (
                "gaps_dist"
                in inspect.signature(smart_gap_threshold_determination).parameters
            ):
                kwargs["gaps_dist"] = msa.site_gappyness
            threshold = smart_gap_threshold_determination(alignment, gaps, **kwargs)
        msa.trim(mode, gap_threshold=threshold)

    elapsed = time.perf_counter() - start
    return {
        "runtime_seconds": elapsed,
        "peak_rss_bytes": _peak_rss_bytes(),
        "output_sha256": _msa_digest(msa, extra),
        "threshold": threshold,
    }


def _execute_worker(case: BenchmarkCase, input_path: Path) -> dict[str, Any]:
    from clipkit.clipkit import execute
    from clipkit.helpers import SeqType
    from clipkit.modes import StopCodonMode, TrimmingMode

    with tempfile.TemporaryDirectory(prefix="clipkit-benchmark-") as temp_dir:
        output_path = Path(temp_dir) / "output"
        auxiliary = (
            str(Path.cwd() / case.auxiliary_file) if case.auxiliary_file else None
        )
        start = time.perf_counter()
        execute(
            input_file=str(input_path),
            input_file_format=case.input_format,
            output_file=str(output_path),
            output_file_format=case.output_format,
            sequence_type=(SeqType(case.sequence_type) if case.sequence_type else None),
            gaps=case.gaps,
            gap_characters=None,
            complement=False,
            codon=case.codon,
            remove_stop_codons=(
                StopCodonMode(case.stop_codons) if case.stop_codons else None
            ),
            ends_only=case.ends_only,
            mode=TrimmingMode[case.mode],
            use_log=False,
            quiet=True,
            dry_run=False,
            validate_only=False,
            report_json=None,
            plot_trim_report=None,
            auxiliary_file=auxiliary,
            threads=case.threads,
        )
        elapsed = time.perf_counter() - start
        output_bytes = output_path.stat().st_size
        output_digest, comparison_type = _output_digest(output_path, case.output_format)

    return {
        "runtime_seconds": elapsed,
        "peak_rss_bytes": _peak_rss_bytes(),
        "output_sha256": output_digest,
        "output_bytes": output_bytes,
        "comparison_type": comparison_type,
    }


def _api_worker(case: BenchmarkCase, input_path: Path) -> dict[str, Any]:
    from clipkit.api import clipkit
    from clipkit.modes import TrimmingMode

    kwargs: dict[str, Any] = {
        "mode": TrimmingMode[case.mode],
        "gaps": case.gaps,
        "sequence_type": case.sequence_type,
        "codon": case.codon,
        "remove_stop_codons": case.stop_codons,
        "ends_only": case.ends_only,
        "threads": case.threads,
    }
    if case.kind == "api_raw":
        kwargs["raw_alignment"] = input_path.read_text()
    else:
        kwargs["input_file_path"] = str(input_path)

    start = time.perf_counter()
    trim_run, stats = clipkit(**kwargs)
    elapsed = time.perf_counter() - start
    return {
        "runtime_seconds": elapsed,
        "peak_rss_bytes": _peak_rss_bytes(),
        "output_sha256": _api_digest(trim_run, stats),
    }


def _cli_worker(case: BenchmarkCase, input_path: Path) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="clipkit-benchmark-") as temp_dir:
        output_path = Path(temp_dir) / "output"
        command = [
            sys.executable,
            "-m",
            "clipkit",
            str(input_path),
            "--output",
            str(output_path),
            "--mode",
            case.mode.replace("_", "-"),
            "--gaps",
            str(case.gaps),
            "--threads",
            str(case.threads),
            "--quiet",
        ]
        if case.sequence_type:
            command.extend(("--sequence_type", case.sequence_type))
        if case.codon:
            command.append("--codon")
        if case.stop_codons:
            command.extend(("--remove_stop_codons", case.stop_codons))
        if case.ends_only:
            command.append("--ends_only")
        if case.input_format:
            command.extend(("--input_file_format", case.input_format))
        if case.output_format:
            command.extend(("--output_file_format", case.output_format))
        if case.auxiliary_file:
            command.extend(("--auxiliary_file", case.auxiliary_file))

        start = time.perf_counter()
        process = subprocess.run(command, capture_output=True, text=True)
        elapsed = time.perf_counter() - start
        if process.returncode != 0:
            raise RuntimeError(
                f"ClipKIT CLI failed ({process.returncode}):\n{process.stderr}"
            )
        output_bytes = output_path.stat().st_size
        output_digest, comparison_type = _output_digest(output_path, case.output_format)

    return {
        "runtime_seconds": elapsed,
        "peak_rss_bytes": _peak_rss_bytes(children=True),
        "output_sha256": output_digest,
        "output_bytes": output_bytes,
        "comparison_type": comparison_type,
    }


def run_worker(case_name: str, input_path: Path) -> int:
    case = CASES[case_name]
    if case.kind == "algorithm":
        result = _algorithm_worker(case, input_path)
    elif case.kind == "execute":
        result = _execute_worker(case, input_path)
    elif case.kind in {"api_path", "api_raw"}:
        result = _api_worker(case, input_path)
    elif case.kind == "cli":
        result = _cli_worker(case, input_path)
    else:
        raise ValueError(f"Unknown benchmark kind: {case.kind}")

    print(json.dumps(result, sort_keys=True))
    return 0


def _invoke_worker(source_root: Path, case: BenchmarkCase, input_path: Path) -> dict:
    worker_command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--case",
        case.name,
        "--input",
        str(input_path),
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [str(source_root), env.get("PYTHONPATH")])
    )
    process = subprocess.run(
        worker_command,
        cwd=source_root,
        env=env,
        capture_output=True,
        text=True,
    )
    if process.returncode != 0:
        raise RuntimeError(
            f"Benchmark worker failed for {case.name} ({process.returncode}):\n"
            f"{process.stderr}"
        )
    return json.loads(process.stdout)


def summarize_samples(case: BenchmarkCase, samples: list[dict]) -> dict:
    digests = {sample["output_sha256"] for sample in samples}
    if len(digests) != 1:
        raise RuntimeError(f"{case.name} output changed between repetitions")

    thresholds = {sample.get("threshold") for sample in samples}
    if len(thresholds) != 1:
        raise RuntimeError(f"{case.name} threshold changed between repetitions")

    runtimes = [sample["runtime_seconds"] for sample in samples]
    peak_rss_values = [
        sample["peak_rss_bytes"]
        for sample in samples
        if sample["peak_rss_bytes"] is not None
    ]
    result = {
        **asdict(case),
        "repetitions": len(samples),
        "runtime_seconds_median": statistics.median(runtimes),
        "runtime_seconds_min": min(runtimes),
        "runtime_seconds_max": max(runtimes),
        "runtime_seconds_range": max(runtimes) - min(runtimes),
        "peak_rss_bytes_median": (
            statistics.median(peak_rss_values) if peak_rss_values else None
        ),
        "peak_rss_bytes_min": min(peak_rss_values) if peak_rss_values else None,
        "peak_rss_bytes_max": max(peak_rss_values) if peak_rss_values else None,
        "output_sha256": samples[0]["output_sha256"],
    }
    for optional_key in ("output_bytes", "threshold", "comparison_type"):
        if optional_key in samples[0]:
            result[optional_key] = samples[0][optional_key]
    return result


def run_case(
    source_root: Path,
    case: BenchmarkCase,
    input_path: Path,
    warmups: int,
    repetitions: int,
) -> dict:
    for _ in range(warmups):
        _invoke_worker(source_root, case, input_path)
    samples = [
        _invoke_worker(source_root, case, input_path) for _ in range(repetitions)
    ]
    return summarize_samples(case, samples)


def _write_generated_fixtures(directory: Path) -> dict[str, Path]:
    rng = random.Random(8675309)
    sparse_path = directory / "generated_sparse_aa.fasta"
    aa_alphabet = "ACDEFGHIKLMNPQRSTVWY"
    with sparse_path.open("w") as handle:
        for row in range(384):
            sequence = "".join(
                "-" if rng.random() < 0.62 else rng.choice(aa_alphabet)
                for _ in range(5000)
            )
            handle.write(f">sparse_{row}\n{sequence}\n")

    stop_path = directory / "generated_stop_nt.fasta"
    sense_codons = ("ATG", "GCT", "CAA", "TTC", "GGA", "ACC", "CGT")
    stop_codons = ("TAA", "TAG", "TGA")
    with stop_path.open("w") as handle:
        for row in range(384):
            codons = [rng.choice(sense_codons) for _ in range(2000)]
            codons[100 + row % 1700] = stop_codons[row % len(stop_codons)]
            codons[-1 - row % 4] = stop_codons[(row + 1) % len(stop_codons)]
            for index in range(250, 2000, 503):
                if (row + index) % 5 == 0:
                    codons[index] = "---"
            handle.write(f">stop_{row}\n{''.join(codons)}\n")

    return {
        "generated_sparse_aa": sparse_path,
        "generated_stop_nt": stop_path,
    }


def _input_paths(source_root: Path, generated: dict[str, Path]) -> dict[str, Path]:
    return {
        name: generated[name] if relative is None else source_root / relative
        for name, relative in FILES.items()
    }


def _git_revision(source_root: Path) -> str | None:
    process = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=source_root,
        capture_output=True,
        text=True,
    )
    return process.stdout.strip() if process.returncode == 0 else None


def _run_source(
    source_root: Path,
    case_names: tuple[str, ...],
    generated: dict[str, Path],
    warmups: int,
    repetitions: int,
) -> dict:
    paths = _input_paths(source_root, generated)
    results = {}
    for name in case_names:
        case = CASES[name]
        result = run_case(source_root, case, paths[case.dataset], warmups, repetitions)
        results[name] = result
        print(
            f"{source_root.name}/{name}: "
            f"{result['runtime_seconds_median']:.6f}s median",
            flush=True,
        )
    _verify_equivalence_groups(results)
    return {
        "source_root": str(source_root),
        "git_revision": _git_revision(source_root),
        "cases": results,
    }


def _verify_equivalence_groups(results: dict[str, dict]) -> None:
    for group in EQUIVALENCE_GROUPS:
        available = [name for name in group if name in results]
        digests = {results[name]["output_sha256"] for name in available}
        if len(digests) > 1:
            raise RuntimeError(
                "Equivalent cases produced different outputs: " + ", ".join(available)
            )


def compare_sources(candidate: dict, reference: dict) -> dict:
    comparison = {}
    for name, candidate_case in candidate["cases"].items():
        reference_case = reference["cases"][name]
        for key in (
            "output_sha256",
            "output_bytes",
            "threshold",
            "comparison_type",
        ):
            if candidate_case.get(key) != reference_case.get(key):
                raise RuntimeError(
                    f"{name} differs from the reference source for {key}: "
                    f"{candidate_case.get(key)!r} != {reference_case.get(key)!r}"
                )
        candidate_time = candidate_case["runtime_seconds_median"]
        reference_time = reference_case["runtime_seconds_median"]
        candidate_memory = candidate_case["peak_rss_bytes_median"]
        reference_memory = reference_case["peak_rss_bytes_median"]
        comparison[name] = {
            "speedup": reference_time / candidate_time,
            "runtime_change_percent": ((candidate_time / reference_time) - 1.0) * 100.0,
            "peak_rss_change_percent": (
                ((candidate_memory / reference_memory) - 1.0) * 100.0
                if candidate_memory is not None and reference_memory
                else None
            ),
        }
    return comparison


def _dependency_versions() -> dict[str, str | None]:
    versions = {}
    for name in ("numpy", "biopython"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="benchmark_smoke.json")
    parser.add_argument("--source-root", type=Path, default=ROOT)
    parser.add_argument(
        "--compare-root",
        type=Path,
        help="reference source tree used for exact-output and performance comparison",
    )
    parser.add_argument(
        "--suite", choices=("smoke", "full", "comprehensive"), default="smoke"
    )
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--case", choices=tuple(CASES), help=argparse.SUPPRESS)
    parser.add_argument("--input", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker:
        if args.case is None or args.input is None:
            parser.error("--worker requires --case and --input")
        return run_worker(args.case, args.input)
    if args.repetitions < 1:
        parser.error("--repetitions must be at least 1")
    if args.warmups < 0:
        parser.error("--warmups cannot be negative")

    case_names = {
        "smoke": SMOKE_CASE_NAMES,
        "full": FULL_CASE_NAMES,
        "comprehensive": COMPREHENSIVE_CASE_NAMES,
    }[args.suite]
    source_root = args.source_root.resolve()
    compare_root = args.compare_root.resolve() if args.compare_root else None

    with tempfile.TemporaryDirectory(prefix="clipkit-benchmark-fixtures-") as temp_dir:
        generated = _write_generated_fixtures(Path(temp_dir))
        candidate = _run_source(
            source_root, case_names, generated, args.warmups, args.repetitions
        )
        reference = (
            _run_source(
                compare_root, case_names, generated, args.warmups, args.repetitions
            )
            if compare_root
            else None
        )

    payload = {
        "timestamp_unix": int(time.time()),
        "suite": args.suite,
        "warmups": args.warmups,
        "repetitions": args.repetitions,
        "python": sys.version,
        "platform": platform.platform(),
        "dependencies": _dependency_versions(),
        "candidate": candidate,
        "reference": reference,
        "candidate_vs_reference": (
            compare_sources(candidate, reference) if reference else None
        ),
    }
    output_path = Path(args.output)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(f"Wrote benchmark report to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

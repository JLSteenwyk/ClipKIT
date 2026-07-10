#!/usr/bin/env python3
import argparse
import hashlib
import inspect
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

FILES = {
    "small": "tests/integration/samples/EOG091N44M8_aa.fa",
    "medium": "tests/integration/samples/EOG092C4VOX_aa_aln.fasta",
    "large": "tests/integration/samples/EOG092C0CZK_aa_aln.fasta",
}

SMOKE_CASES = [
    ("gappy_small", "end_to_end", "small", "gappy"),
    ("smart_gap_medium", "end_to_end", "medium", "smart_gap"),
    ("kpic_large", "end_to_end", "large", "kpic"),
]

FULL_CASES = [
    ("construct_small", "algorithm", "small", "construct"),
    ("construct_medium", "algorithm", "medium", "construct"),
    ("construct_large", "algorithm", "large", "construct"),
    ("gappy_large_algorithm", "algorithm", "large", "gappy"),
    ("smart_gap_large_algorithm", "algorithm", "large", "smart_gap"),
    ("kpic_large_algorithm", "algorithm", "large", "kpic"),
    ("entropy_large_algorithm", "algorithm", "large", "entropy"),
    (
        "composition_bias_large_algorithm",
        "algorithm",
        "large",
        "composition_bias",
    ),
    ("gappy_small", "end_to_end", "small", "gappy"),
    ("smart_gap_medium", "end_to_end", "medium", "smart_gap"),
    ("gappy_large", "end_to_end", "large", "gappy"),
    ("kpic_large", "end_to_end", "large", "kpic"),
    ("c3_large", "end_to_end", "large", "c3"),
]


def _peak_rss_bytes() -> int | None:
    try:
        import resource
    except ImportError:
        return None

    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak_rss if sys.platform == "darwin" else peak_rss * 1024


def _msa_digest(msa) -> str:
    payload = (
        msa.seq_records.tobytes()
        + msa._site_positions_to_keep.astype("int64").tobytes()
        + msa._site_positions_to_trim.astype("int64").tobytes()
    )
    return hashlib.sha256(payload).hexdigest()


def run_worker(kind: str, operation: str, input_path: Path) -> int:
    from clipkit.modes import TrimmingMode

    if kind == "algorithm":
        from Bio import AlignIO

        from clipkit.msa import MSA
        from clipkit.settings import DEFAULT_AA_GAP_CHARS
        from clipkit.smart_gap_helper import smart_gap_threshold_determination

        alignment = AlignIO.read(input_path, "fasta")
        start = time.perf_counter()
        msa = MSA.from_bio_msa(alignment, DEFAULT_AA_GAP_CHARS)
        threshold = None

        if operation == "gappy":
            msa.trim(TrimmingMode.gappy, gap_threshold=0.9)
        elif operation == "smart_gap":
            kwargs = {"seq_records": msa.seq_records}
            if "gaps_dist" in inspect.signature(
                smart_gap_threshold_determination
            ).parameters:
                kwargs["gaps_dist"] = msa.site_gappyness
            threshold = smart_gap_threshold_determination(
                alignment,
                DEFAULT_AA_GAP_CHARS,
                **kwargs,
            )
            msa.trim(TrimmingMode.smart_gap, gap_threshold=threshold)
        elif operation == "kpic":
            msa.trim(TrimmingMode.kpic, gap_threshold=0.9)
        elif operation == "entropy":
            msa.trim(TrimmingMode.entropy, gap_threshold=0.8)
        elif operation == "composition_bias":
            msa.trim(TrimmingMode.composition_bias, gap_threshold=0.8)

        elapsed = time.perf_counter() - start
        result = {
            "runtime_seconds": elapsed,
            "peak_rss_bytes": _peak_rss_bytes(),
            "output_sha256": _msa_digest(msa),
            "threshold": threshold,
        }
    else:
        from clipkit.clipkit import execute

        with tempfile.TemporaryDirectory(prefix="clipkit-benchmark-") as temp_dir:
            output_path = Path(temp_dir) / "output.fa"
            start = time.perf_counter()
            execute(
                input_file=str(input_path),
                input_file_format=None,
                output_file=str(output_path),
                output_file_format=None,
                sequence_type=None,
                gaps=0.9,
                gap_characters=None,
                complement=False,
                codon=False,
                ends_only=False,
                mode=getattr(TrimmingMode, operation),
                use_log=False,
                quiet=True,
                dry_run=False,
                validate_only=False,
                report_json=None,
                plot_trim_report=None,
                auxiliary_file=None,
                threads=1,
            )
            elapsed = time.perf_counter() - start
            output = output_path.read_bytes()

        result = {
            "runtime_seconds": elapsed,
            "peak_rss_bytes": _peak_rss_bytes(),
            "output_sha256": hashlib.sha256(output).hexdigest(),
            "output_bytes": len(output),
        }

    print(json.dumps(result))
    return 0


def run_case(
    source_root: Path,
    kind: str,
    scale: str,
    operation: str,
    repetitions: int,
) -> dict:
    input_path = source_root / FILES[scale]
    worker_command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--kind",
        kind,
        "--operation",
        operation,
        "--input",
        str(input_path),
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        filter(None, [str(source_root), env.get("PYTHONPATH")])
    )

    samples = []
    for _ in range(repetitions):
        proc = subprocess.run(
            worker_command,
            cwd=source_root,
            env=env,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"Benchmark worker failed ({proc.returncode}):\n{proc.stderr}"
            )
        samples.append(json.loads(proc.stdout))

    digests = {sample["output_sha256"] for sample in samples}
    if len(digests) != 1:
        raise RuntimeError("Benchmark output changed between repetitions")

    runtimes = [sample["runtime_seconds"] for sample in samples]
    peak_rss_values = [
        sample["peak_rss_bytes"]
        for sample in samples
        if sample["peak_rss_bytes"] is not None
    ]
    result = {
        "kind": kind,
        "scale": scale,
        "operation": operation,
        "repetitions": repetitions,
        "runtime_seconds_median": statistics.median(runtimes),
        "runtime_seconds_min": min(runtimes),
        "runtime_seconds_max": max(runtimes),
        "peak_rss_bytes_median": (
            statistics.median(peak_rss_values) if peak_rss_values else None
        ),
        "output_sha256": samples[0]["output_sha256"],
    }
    for optional_key in ("output_bytes", "threshold"):
        if optional_key in samples[0]:
            result[optional_key] = samples[0][optional_key]
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="benchmark_smoke.json")
    parser.add_argument("--source-root", type=Path, default=ROOT)
    parser.add_argument("--suite", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--kind", choices=("algorithm", "end_to_end"))
    parser.add_argument("--operation")
    parser.add_argument("--input", type=Path)
    args = parser.parse_args()

    if args.worker:
        return run_worker(args.kind, args.operation, args.input)
    if args.repetitions < 1:
        parser.error("--repetitions must be at least 1")

    source_root = args.source_root.resolve()
    cases = SMOKE_CASES if args.suite == "smoke" else FULL_CASES
    results = {}
    for name, kind, scale, operation in cases:
        results[name] = run_case(
            source_root,
            kind,
            scale,
            operation,
            args.repetitions,
        )
        print(
            f"{name}: {results[name]['runtime_seconds_median']:.6f}s median",
            flush=True,
        )

    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=source_root,
        capture_output=True,
        text=True,
    )
    payload = {
        "timestamp_unix": int(time.time()),
        "source_root": str(source_root),
        "git_revision": revision.stdout.strip() if revision.returncode == 0 else None,
        "suite": args.suite,
        "python": sys.version,
        "platform": platform.platform(),
        "cases": results,
    }
    output_path = Path(args.output)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(f"Wrote benchmark report to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

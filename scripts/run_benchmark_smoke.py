#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_case(name: str, cmd: list[str]) -> dict:
    start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = round(time.perf_counter() - start, 6)
    return {
        "name": name,
        "command": " ".join(cmd),
        "returncode": proc.returncode,
        "runtime_seconds": elapsed,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="benchmark_smoke.json",
        help="Path to write benchmark JSON output.",
    )
    args = parser.parse_args()

    sample = ROOT / "tests" / "integration" / "samples" / "EOG091N44M8_aa.fa"
    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)

    python = sys.executable
    base = [python, "-m", "clipkit-runner", str(sample), "-m", "gappy"]

    cases = [
        ("gappy_threads_1", base + ["-t", "1", "-o", str(out_dir / "bench.gappy.t1.fa")]),
        ("gappy_threads_4", base + ["-t", "4", "-o", str(out_dir / "bench.gappy.t4.fa")]),
    ]

    results = [run_case(name, cmd) for name, cmd in cases]
    payload = {
        "timestamp_unix": int(time.time()),
        "cases": results,
    }

    output_path = Path(args.output)
    output_path.write_text(json.dumps(payload, indent=2))

    failed = [case for case in results if case["returncode"] != 0]
    if failed:
        print("Benchmark smoke run failed:")
        for case in failed:
            print(f"- {case['name']} returncode={case['returncode']}")
        return 1

    print(f"Wrote benchmark report to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

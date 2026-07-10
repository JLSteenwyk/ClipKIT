# ClipKIT performance optimization report

## Scope

Performance work was conducted on `codex/faster-msa-construction` and compared
with `master` at `e0d93b4`. The optimized code benchmarked below ends at
`52814f0`.

The implementation preserves ClipKIT's public APIs, trimming decisions,
sequence ordering, headers, thresholds, and output formatting. No trimming
algorithm or parameter was changed.

## Optimizations

1. **Fast MSA construction (`47319c4`)**
   converts complete sequence strings directly into a NumPy Unicode matrix,
   avoiding a Python object for every input character.
2. **Batched column statistics (`ad1314a`)**
   counts compact character codes across batches of columns instead of sorting
   every column independently. Counts use the smallest safe unsigned dtype and
   retain a compact-Unicode fallback for unusually wide alphabets.
3. **Shared smart-gap statistics (`4dc2046`)**
   feeds the cached per-site gap distribution into smart-gap threshold
   selection, removing a duplicate full-alignment scan.
4. **Fast output materialization (`52814f0`)**
   views contiguous NumPy rows as complete Unicode strings instead of creating
   millions of one-character Python objects before writing.

## Benchmark method

Benchmarks ran on `macOS-26.4.1-arm64-arm-64bit` with Python 3.11, NumPy
1.26.4, and Biopython 1.83. Each result is the median of fresh processes: seven
repetitions for the small case and three for medium and large cases. Peak RSS is
the median process maximum.

The representative alignments were:

| Scale | Sequences | Sites | Alignment cells |
| --- | ---: | ---: | ---: |
| Small | 117 | 1,026 | 120,042 |
| Medium | 1,480 | 12,977 | 19,205,960 |
| Large | 1,478 | 29,838 | 44,100,564 |

Algorithm benchmarks load the Biopython alignment before timing, then include
MSA construction and the named trimming calculation. Their equality check
hashes the sequence matrix and exact keep/trim position arrays.

| Algorithm case | Baseline | Optimized | Speedup | Peak RSS baseline → optimized |
| --- | ---: | ---: | ---: | ---: |
| Construct, small | 0.011783 s | 0.000221 s | 53.31× | 57.4 → 56.2 MiB |
| Construct, medium | 1.695093 s | 0.019756 s | 85.80× | 456.2 → 313.7 MiB |
| Construct, large | 3.826577 s | 0.047105 s | 81.23× | 950.3 → 631.6 MiB |
| Gappy, large | 4.047438 s | 0.192978 s | 20.97× | 1,034.3 → 661.2 MiB |
| Smart-gap, large | 4.319067 s | 0.194392 s | 22.22× | 1,034.9 → 662.0 MiB |
| KPIC, large | 4.636736 s | 0.193192 s | 24.00× | 975.4 → 663.2 MiB |
| Entropy, large | 4.620146 s | 0.383034 s | 12.06× | 951.1 → 671.2 MiB |
| Composition-bias, large | 4.680238 s | 0.327606 s | 14.29× | 975.6 → 662.4 MiB |

End-to-end benchmarks include FASTA parsing, trimming, materialization, and
FASTA writing. Equality is a byte-for-byte SHA-256 comparison of output files.

| End-to-end case | Baseline | Optimized | Speedup | Peak RSS baseline → optimized |
| --- | ---: | ---: | ---: | ---: |
| Gappy, small | 0.022275 s | 0.004334 s | 5.14× | 57.6 → 55.6 MiB |
| Smart-gap, medium | 2.555377 s | 0.404308 s | 6.32× | 530.5 → 300.3 MiB |
| Gappy, large | 4.303231 s | 0.402573 s | 10.69× | 733.5 → 333.7 MiB |
| KPIC, large | 5.194196 s | 0.515740 s | 10.07× | 735.2 → 414.6 MiB |
| C3, large | 4.794962 s | 0.490340 s | 9.78× | 837.7 → 438.8 MiB |

Every comparison above produced identical hashes. No benchmarked workload was
slower after optimization.

## Correctness verification

- The full suite passes: 345 tests on the supported NumPy 1.26.4/Biopython
  1.83 combination.
- New regression tests compare batched counts, frequencies, gappyness,
  entropy, composition bias, and classifications with the original
  per-column calculations, including mixed case, custom gaps, empty rows, and
  wide Unicode alphabets.
- Thirteen non-auxiliary CLI modes produced byte-identical FASTA output between
  baseline and optimized revisions on the same input.
- Smart-gap thresholds were exactly equal in all benchmark comparisons.
- Existing integration tests continue to verify supported output formats,
  complementary output, logging, eComp behavior, headers, and ordering.

## Reproducing benchmarks

The scheduled smoke benchmark now covers small, medium, and large workloads,
records output hashes and peak RSS, and can benchmark another worktree:

```shell
git worktree add /tmp/clipkit-baseline e0d93b4
python scripts/run_benchmark_smoke.py \
  --source-root /tmp/clipkit-baseline \
  --suite full \
  --output benchmark-baseline.json
python scripts/run_benchmark_smoke.py \
  --source-root . \
  --suite full \
  --output benchmark-optimized.json
```

Runtime varies by machine, so compare reports produced on the same otherwise
idle system and interpreter.

## Remaining bottlenecks

- Biopython FASTA parsing and format-specific writing now dominate common
  large end-to-end runs. Replacing them would put support for all current
  formats and exact formatting at risk, so this work leaves them intact.
- Heterotachy mode is dominated by Biopython's parsimony guide-tree search. On
  the 24-sequence, 231-site fixture, tree construction consumed 7.21 of 7.31
  profiled seconds. Substituting a faster tree search could change the guide
  tree and trimming output, so it was not done.
- Entropy's exact per-site floating-point calculation remains more expensive
  than count-based classification. Vectorizing it could alter values at
  four-decimal rounding boundaries and has a poor correctness-to-reward ratio
  now that the surrounding count path is substantially faster.

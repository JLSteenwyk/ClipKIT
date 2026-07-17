# ClipKIT performance optimization report

## Scope

The current performance work compares `master` at `fb95e1b` with the runtime
optimization series ending at `701009b`. The benchmark and differential-test
commits are `2afba34` and `0610ca5`.

The implementation preserves ClipKIT's public APIs, trimming decisions,
sequence ordering, headers, thresholds, statistics, reports, logs, warnings,
errors, and output formatting. No trimming algorithm, biological rule, or
parameter was changed, and no runtime dependency was added.

## Current optimizations

1. **Stop-codon case handling (`aa9e9ca`)** uses the MSA's existing mixed-case
   metadata to avoid applying Unicode uppercase conversion to every matrix cell
   when the input is already uppercase. An independent working array is retained
   so allocation and mutation semantics remain unchanged.
2. **Vectorized site scoring (`007a2f0`)** evaluates entropy and composition-bias
   over all columns for each small alphabet state instead of looping over every
   alignment column in Python. Composition-bias columns are grouped by observed
   state count to preserve NumPy's original reduction order at four-decimal
   rounding boundaries.
3. **Vectorized codon expansion (`701009b`)** expands unique codon blocks in one
   NumPy operation instead of constructing every triplet in a Python loop.

These changes build on the MSA construction, batched column counting, shared
smart-gap statistics, and output materialization improvements released in
ClipKIT 2.12.1 and 2.12.2.

## Benchmark method

Benchmarks ran on Apple M2 macOS 26.4.1 with Python 3.11.14, NumPy 1.26.4,
and Biopython 1.83. Every case used one discarded warm-up followed by five
fresh-process measurements. Reports include median/range wall-clock time,
CPU time, and peak RSS. Tables below show wall time. Algorithm equality hashes
the sequence matrix plus
the exact keep/trim arrays; end-to-end deterministic formats use byte-for-byte
SHA-256. ECOMP uses a decoded alignment/metadata digest because its gzip
fallback embeds a creation timestamp.

The comprehensive suite contains 50 cases spanning small, medium, and large
amino-acid and nucleotide inputs, deterministic sparse/gappy inputs, every
trimming mode, codon trimming, all three stop-codon modes, CLI and API entry
points, one and four requested threads, ECOMP input, and all supported output
formats.

### Optimized hot paths

| Case | Baseline median (range) | Optimized median (range) | Speedup |
| --- | ---: | ---: | ---: |
| Entropy, large algorithm | 0.3375 s (0.3343–0.3461) | 0.1833 s (0.1755–0.1951) | 1.84× |
| Composition-bias, large algorithm | 0.2864 s (0.2832–0.2924) | 0.1751 s (0.1725–0.1989) | 1.64× |
| Terminal stop masking, algorithm | 0.3817 s (0.3691–0.3917) | 0.0499 s (0.0339–0.0548) | 7.64× |
| Internal stop masking, algorithm | 0.3793 s (0.3685–0.3892) | 0.0580 s (0.0416–0.0788) | 6.54× |
| All stop masking, algorithm | 0.3738 s (0.3646–0.3875) | 0.0524 s (0.0412–0.0772) | 7.13× |
| Terminal stop masking, end to end | 0.4116 s (0.3947–0.4210) | 0.0882 s (0.0834–0.1023) | 4.67× |
| Internal stop masking, end to end | 0.3982 s (0.3978–0.4174) | 0.0852 s (0.0762–0.0917) | 4.67× |
| All stop masking, end to end | 0.4071 s (0.3928–0.4145) | 0.0948 s (0.0775–0.1172) | 4.29× |

Unchanged large common paths stayed within noise: gappy end-to-end was
0.3737 s before and 0.3682 s after; KPIC end-to-end was 0.4990 s before and
0.4840 s after. A separate 12-repetition interleaved audit of cases that looked
slower in sequential reports found no regression: large KPIC algorithm −2.2%,
medium smart-gap −0.5%, sparse gappy −4.8%, and small gappy −14.7%.

### Memory

Large entropy and composition-bias median process peak RSS changed by +2.7%,
inside the observed run-to-run range; targeted interleaved measurements ranged
from −2.4% to −0.2%. Other unchanged large workflows stayed within about 1%.

macOS `ru_maxrss` reports 12–17% higher residency for the stop-codon
end-to-end cases after they became more than four times faster. Allocation
tracing repeated five times measured the same 28,162,648-byte peak before and
after, and the optimized masking function itself allocates no more memory than
the baseline. This is a residency/sampling effect from the much shorter-lived
process, not a memory-for-speed tradeoff; it is disclosed here because the
benchmark report intentionally retains raw process RSS.

## Correctness verification

- The comprehensive candidate and baseline reports have identical hashes,
  output byte lengths, inferred thresholds, and comparison modes for all 50
  cases.
- Randomized scalar-oracle tests cover entropy, composition bias, and every
  stop-codon mode, including lowercase input, trailing gap codons, and mixed
  complete/incomplete codons.
- Codon-expansion tests cover partial final codons, duplicates, unsorted site
  input, and minimal alignments.
- Existing integration tests continue to cover supported formats, CLI/API
  behavior, complementary output, logs, reports, headers, ordering, cached
  statistic invalidation, and stop-codon summaries.

## Reproducing benchmarks

Use a clean baseline worktree and the supported dependency versions:

```shell
git worktree add --detach /tmp/clipkit-baseline fb95e1b
python scripts/run_benchmark_smoke.py \
  --source-root . \
  --compare-root /tmp/clipkit-baseline \
  --suite comprehensive \
  --warmups 1 \
  --repetitions 5 \
  --output benchmark-comparison.json
```

The runner fails if retained repetitions differ or if the candidate and
reference outputs, thresholds, or lengths differ. Runtime varies by machine,
so compare trees using the same idle system and interpreter.

## Remaining bottlenecks

- Biopython FASTA parsing and format-specific writing dominate common large
  end-to-end runs. Replacing them risks changing supported-format behavior and
  exact serialization, so they remain intact.
- Batched character counting dominates common gappy and KPI/KPIC algorithm
  paths. Batch sizes from 128 through 4,096 were measured; the existing 1,024
  setting remained fastest on the large fixture.
- Heterotachy remains dominated by Biopython parsimony guide-tree search.
  Substituting another tree algorithm could change the guide tree and trimming
  result, so it was not attempted.

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
| Entropy, large algorithm | 0.3397 s (0.3348–0.3550) | 0.1780 s (0.1727–0.2335) | 1.91× |
| Composition-bias, large algorithm | 0.2865 s (0.2849–0.2916) | 0.1751 s (0.1668–0.1818) | 1.64× |
| Terminal stop masking, algorithm | 0.3757 s (0.3714–0.3806) | 0.0543 s (0.0455–0.0695) | 6.92× |
| Internal stop masking, algorithm | 0.3852 s (0.3684–0.3878) | 0.0588 s (0.0380–0.0697) | 6.55× |
| All stop masking, algorithm | 0.3639 s (0.3551–0.4021) | 0.0524 s (0.0384–0.0565) | 6.94× |
| Terminal stop masking, end to end | 0.4042 s (0.4003–0.4145) | 0.0837 s (0.0799–0.1027) | 4.83× |
| Internal stop masking, end to end | 0.4049 s (0.4002–0.4128) | 0.0922 s (0.0794–0.1064) | 4.39× |
| All stop masking, end to end | 0.3990 s (0.3950–0.4144) | 0.0975 s (0.0903–0.1161) | 4.09× |

Unchanged large common paths stayed within noise: gappy end-to-end was
0.3768 s before and 0.3672 s after; KPIC end-to-end was 0.4738 s before and
0.4760 s after. A separate 12-repetition interleaved audit of cases that looked
slower in sequential reports found no regression: large KPIC algorithm −2.2%,
medium smart-gap −0.5%, sparse gappy −4.8%, and small gappy −14.7%.

### Memory

Large entropy and composition-bias median process peak RSS changed by −0.2%
and +0.2%. Other unchanged large workflows stayed within about 2%.

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

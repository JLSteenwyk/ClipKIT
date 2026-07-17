# ClipKIT performance optimization report

## Scope

The current performance work has two measured baselines. The earlier round
compares `master` at `fb95e1b` with the runtime series ending at `701009b`.
The focused core-mode round compares `9f456cb` with the series ending at
`f3b0c02`; its benchmark and differential-test commits are `3b8c571` and
`e2368ba`.

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
4. **Combined count offsets (`aa95fbd`)** folds the character-code origin into
   the column offsets, removing one full in-place pass over every counting
   batch.
5. **Shared combined-mode counts (`a1f2b4e`)** calculates classification before
   gappyness for uppercase KPI/KPIC combined modes and derives both statistics
   from the same exact count matrix. Mixed-case alignments retain separate raw
   and normalized passes, preserving case and memory behavior.
6. **Fused count encoding (`f3b0c02`)** casts Unicode code points and applies
   column offsets in one NumPy operation instead of materializing and then
   revisiting the temporary integer matrix.

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

The comprehensive suite contains 104 cases spanning small, medium, and large
amino-acid and nucleotide inputs, deterministic sparse/gappy inputs, every
trimming mode, codon trimming, all three stop-codon modes, CLI and API entry
points, one and four requested threads, ECOMP input, and all supported output
formats.

The focused `core` suite contains 54 cases. It exercises gappy, smart-gap, KPI,
KPIC, and all four gap/classification combinations on a large sparse AA input,
a medium sparse NT input, and deterministic dense AA and NT inputs. It also
covers small and large end-to-end runs, CLI/API entry points, and requested
thread counts of one and four. Algorithm/API results independently hash exact
keep/trim arrays, classification arrays, and gap distributions and record
resolved thresholds and effective threads.

### Focused core-mode results

The tables compare the final implementation with `9f456cb`. The large AA
alignment contains 1,478 sequences, 29,838 columns, and about 95% gap
characters. Ranges are the minimum and maximum of five retained fresh-process
measurements.

| Algorithm case | Baseline median (range) | Optimized median (range) | Speedup | Median RSS change |
| --- | ---: | ---: | ---: | ---: |
| gappy | 0.2092 s (0.1909–0.4147) | 0.1658 s (0.1561–0.1827) | 1.26× | −2.6% |
| smart-gap | 0.1840 s (0.1781–0.2122) | 0.1656 s (0.1572–0.1803) | 1.11× | −2.4% |
| KPI | 0.1882 s (0.1751–0.2305) | 0.1585 s (0.1523–0.1705) | 1.19× | +2.1% |
| KPIC | 0.1862 s (0.1738–0.2010) | 0.1645 s (0.1571–0.2027) | 1.13× | +2.3% |
| KPI-gappy | 0.3806 s (0.3290–0.3958) | 0.1679 s (0.1600–0.1808) | 2.27× | −3.6% |
| KPIC-gappy | 0.3426 s (0.3302–0.4592) | 0.1598 s (0.1577–0.1929) | 2.14× | −3.8% |
| KPI-smart-gap | 0.3446 s (0.3333–0.3524) | 0.1605 s (0.1558–0.1689) | 2.15× | −4.1% |
| KPIC-smart-gap | 0.3296 s (0.3181–0.3404) | 0.1687 s (0.1547–0.1986) | 1.95× | −1.3% |

| End-to-end case | Baseline median (range) | Optimized median (range) | Speedup | Median RSS change |
| --- | ---: | ---: | ---: | ---: |
| gappy | 0.3377 s (0.3286–0.3843) | 0.3144 s (0.3086–0.3204) | 1.07× | +1.2% |
| smart-gap | 0.6198 s (0.6166–0.6366) | 0.6138 s (0.6018–0.6345) | 1.01× | −0.4% |
| KPI | 0.4134 s (0.4120–0.4192) | 0.4035 s (0.4002–0.4144) | 1.02× | +0.1% |
| KPIC | 0.4534 s (0.4455–0.4598) | 0.4228 s (0.4158–0.4304) | 1.07× | −0.0% |
| KPI-gappy | 0.4728 s (0.4640–0.4951) | 0.3155 s (0.3141–0.3243) | 1.50× | −3.2% |
| KPIC-gappy | 0.4647 s (0.4587–0.4896) | 0.3161 s (0.3157–0.3175) | 1.47× | −1.5% |
| KPI-smart-gap | 0.5526 s (0.5368–0.5820) | 0.3956 s (0.3941–0.4091) | 1.40× | −1.2% |
| KPIC-smart-gap | 0.5747 s (0.5661–0.5807) | 0.4227 s (0.4190–0.4298) | 1.36× | +0.4% |

The final sequential report made several 3–90 ms cases appear slower because
the two source trees run in separate blocks. A 12-repetition alternating audit
reversed all material apparent regressions: API raw KPIC-smart-gap was 1.23×
faster, dense-NT smart-gap 1.04×, medium-NT KPIC 1.07×, dense-NT KPIC 1.30×,
CLI smart-gap 1.01×, and dense-AA KPI 1.19×. A subsequent 30-repetition audit
of medium-NT KPI measured a 1.3% median and 3.2% mean wall-time improvement.

### Profiling and optimization decisions

At `9f456cb`, deterministic `cProfile` runs attributed 0.136–0.168 seconds of
each large base-mode call to `_column_character_counts`. Combined modes called
that function twice (0.259–0.269 seconds total), while MSA construction used
about 0.04 seconds and smart-gap threshold selection about 0.001 seconds. This
made count sharing the only material combined-mode opportunity.

Each accepted change was benchmarked against the immediately preceding commit:

- Combining offset passes improved the large paths by 1.03–1.10× CPU time.
- Sharing classification/gap counts improved large combined modes by
  1.67–1.77× without changing base modes.
- Fusing cast and offset improved all large modes by 1.04–1.12×. Twenty
  interleaved kernel repetitions measured 1.05× on dense AA and 1.54× on dense
  NT. Allocation tracing showed only a 132 KB transient NumPy iterator buffer,
  0.49% of the 27 MB counting peak, with no retained allocation.

Rejected or deferred ideas included 32-bit encoded temporaries (not faster in
interleaved tests), another batch-size change (1,024 remained best in the prior
128–4,096 sweep), vectorizing sub-millisecond smart-gap slope/mask operations,
sharing raw counts for mixed-case inputs (case normalization and peak-memory
risk), and replacing Biopython parsing/writing (high compatibility risk).

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

Focused large combined algorithm cases reduced median RSS by 1.3–4.1%. Large
base-mode median RSS varied from −2.6% to +2.3%, within process-residency noise;
large end-to-end cases varied from −0.4% to +1.2% except for the combined modes,
which were mostly lower. Direct allocation tracing is described above.

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

- The focused candidate and `9f456cb` reports have identical biological output
  hashes, byte lengths, thresholds, exact keep/trim arrays, classification
  arrays, gap distributions, and effective thread counts for all 54 cases.
- Randomized scalar-oracle tests cover all eight core modes across mixed-case
  inputs and exact `>`/`>=` threshold boundaries. Dedicated tests prove that
  uppercase combined modes use one count pass and mixed-case modes retain the
  separate raw/normalized passes.
- The final comprehensive candidate and baseline reports have identical hashes,
  output byte lengths, inferred thresholds, exact statistic/position metadata,
  comparison modes, and thread decisions for all 104 cases.
- Randomized scalar-oracle tests cover entropy, composition bias, and every
  stop-codon mode, including lowercase input, trailing gap codons, and mixed
  complete/incomplete codons.
- Codon-expansion tests cover partial final codons, duplicates, unsorted site
  input, and minimal alignments.
- Existing integration tests continue to cover supported formats, CLI/API
  behavior, complementary output, logs, reports, headers, ordering, cached
  statistic invalidation, and stop-codon summaries.

## Reproducing benchmarks

Use a clean baseline worktree and the supported dependency versions. For the
focused core-mode comparison:

```shell
git worktree add --detach /tmp/clipkit-core-baseline 9f456cb
python scripts/run_benchmark_smoke.py \
  --source-root . \
  --compare-root /tmp/clipkit-core-baseline \
  --suite core \
  --warmups 1 \
  --repetitions 5 \
  --output core-modes-final-vs-9f456cb.json
```

The runner fails if retained repetitions differ or if the candidate and
reference outputs, thresholds, arrays, thread decisions, or lengths differ.
Runtime varies by machine, so compare trees using the same idle system and
interpreter. The retained comprehensive/focused JSON reports, compact focused
CSV, and interleaved audit are in `benchmark-results/`.

## Remaining bottlenecks

- Biopython FASTA parsing and format-specific writing dominate common large
  end-to-end runs. Replacing them risks changing supported-format behavior and
  exact serialization, so they remain intact.
- Batched character counting still dominates common gappy and KPI/KPIC
  algorithm paths after the safe encoding improvements. Further changes would
  need to outperform the current 1,024-column kernel across sparse/dense AA and
  NT shapes without increasing its temporary peak.
- Heterotachy remains dominated by Biopython parsimony guide-tree search.
  Substituting another tree algorithm could change the guide tree and trimming
  result, so it was not attempted.

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
The remaining unchanged non-core outlier, internal stop masking, was 3.3%
faster by median and 0.7% faster by mean in its own 30-repetition audit.

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

## 1.75× acceleration project

This round compares ClipKIT 2.15.0 at `6044b1e` with the performance branch.
The [targeted samples](benchmark-results/clipkit-175x-targeted.json) identify
the revision pair for each optimization, and the
[interleaved regression audit](benchmark-results/clipkit-175x-regression-audit.json)
compares `161eb0c` directly with the starting revision. Targeted comparisons
use the reference shown in the data. Treat each speedup as an incremental
result against its stated reference.

| Workload | Reference median | Candidate median | Speedup |
| --- | ---: | ---: | ---: |
| CLI version startup | 246 ms | 224 ms | 1.10× |
| Medium gappy CLI with logging | 5.13 s | 1.48 s | 3.47× |
| Medium mixed-case KPIC smart-gap CLI | 4.33 s | 0.58 s | 7.47× |
| Medium heterotachy CLI | 1.07 s | 0.63 s | 1.70× |
| Large gappy end to end, bounded FASTA reader | 538 ms | 504 ms | 1.07× |
| Medium smart-gap end to end, batched FASTA writer | 393 ms | 347 ms | 1.13× |
| Large HTML report CLI | 1.48 s | 1.01 s | 1.47× |
| Large entropy CLI | 1.20 s | 1.06 s | 1.13× |
| Medium ECOMP output | 347 ms | 339 ms | 1.02× |

The regression audit found a specific import tradeoff. Tiny heterotachy's
first in-process run moved from 6.2 ms to 27.4 ms because the timed call now
includes the deferred tree import. The full CLI for the same input moved
from 244 ms to 250 ms. Several tiny output-format cases incurred
0.10–0.23 ms higher medians. In the same audit, large KPIC smart-gap output
improved from 765 ms to 711 ms. The retained changes aim to reduce actual
CLI and API work. A separate prototype for lazy output-statistics formatting
gave mixed results on repeated tiny cases and was discarded.

`make test` passed 520 unit and 135 integration tests on Linux for each of
CPython 3.12.14, 3.13.15, 3.14.7, and 3.15.2 at `161eb0c`. The same suite
also passed on 3.12.14 with the NumPy 1.26.0 and Biopython 1.82 dependency
floors. The CI version check passed on the pushed branch. The full hosted
macOS matrix has not run for this branch.

Run the suites under a single Python 3.12 environment with the benchmark
dependencies recorded in the JSON:

```bash
git worktree add --detach /tmp/clipkit-175x-reference 6044b1e
git switch perf/clipkit-175x
python3.12 -m venv /tmp/clipkit-175x-venv
source /tmp/clipkit-175x-venv/bin/activate
python -m pip install -r requirements.txt . numpy==2.0.2 biopython==1.83
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python scripts/run_benchmark_smoke.py --compare-root /tmp/clipkit-175x-reference \
  --suite core --warmups 1 --repetitions 5 --output /tmp/clipkit-core.json
python scripts/run_benchmark_smoke.py --compare-root /tmp/clipkit-175x-reference \
  --suite comprehensive --warmups 1 --repetitions 5 --output /tmp/clipkit-comprehensive.json
make test
```

The runner compares deterministic bytes and alignment statistics and decodes
ECOMP for its content and metadata check. The comprehensive suite adds CLI
startup, HTML reporting, logs and complementary output, ambiguity-handling
policies, large entropy and composition scoring, and medium heterotachy
and mixed-case fixtures to the earlier coverage.

### Measured checkpoint at `161eb0c`

The [comprehensive JSON](benchmark-results/clipkit-175x-comprehensive-161eb0c.json)
contains fresh-process wall and CPU samples, peak RSS, and the comparison
metadata for every case. Both revisions used CPython 3.12.14, NumPy 2.0.2,
and Biopython 1.83 on Linux x86-64. Five recorded runs followed one warm-up
for each case and each source tree.

### Same-machine suite comparison

The comprehensive suite compared 119 cases exactly against `6044b1e`.
Each row shows a representative case for one public trimming mode.

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| `gappy` | 623.25 ms | 465.41 ms | 1.34× |
| `smart_gap` | 1,013.93 ms | 769.70 ms | 1.32× |
| `block_gappy` | 646.82 ms | 476.31 ms | 1.36× |
| `gappyout` | 1,023.80 ms | 777.55 ms | 1.32× |
| `kpi` | 755.24 ms | 685.23 ms | 1.10× |
| `kpic` | 784.87 ms | 681.81 ms | 1.15× |
| `kpi_gappy` | 646.98 ms | 585.70 ms | 1.10× |
| `kpic_gappy` | 631.18 ms | 587.11 ms | 1.08× |
| `kpi_smart_gap` | 738.11 ms | 664.50 ms | 1.11× |
| `kpic_smart_gap` | 756.53 ms | 695.92 ms | 1.09× |
| `entropy` | 1,005.06 ms | 740.19 ms | 1.36× |
| `composition_bias` | 1,035.69 ms | 752.35 ms | 1.38× |
| `heterotachy` | 1,084.93 ms | 639.86 ms | 1.70× |
| `cst` | 0.69 ms | 0.60 ms | 1.15× |
| `c3` | 6.35 ms | 6.38 ms | 0.99× |

### Major workflows

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| CLI help startup | 237.00 ms | 213.56 ms | 1.11× |
| CLI version startup | 233.53 ms | 224.67 ms | 1.04× |
| CLI small FASTA trimming | 251.82 ms | 222.51 ms | 1.13× |
| Python API from file | 4.38 ms | 4.45 ms | 0.99× |
| Python API from string | 4.71 ms | 4.96 ms | 0.95× |
| HTML report | 1,820.13 ms | 945.85 ms | 1.92× |
| Log and complementary output | 1,190.12 ms | 1,067.97 ms | 1.11× |
| Mixed-case amino-acid CLI | 1,050.96 ms | 275.36 ms | 3.82× |
| ECOMP input | 2.19 ms | 2.22 ms | 0.98× |
| ECOMP output | 403.95 ms | 336.42 ms | 1.20× |
| Codon trimming | 5.22 ms | 5.19 ms | 1.01× |
| Terminal stop-codon masking | 185.19 ms | 183.74 ms | 1.01× |
| Internal stop-codon masking | 181.15 ms | 183.70 ms | 0.99× |
| All stop-codon masking | 180.90 ms | 182.03 ms | 0.99× |
| One requested thread | 463.95 ms | 346.58 ms | 1.34× |
| Four requested threads | 445.84 ms | 348.71 ms | 1.28× |
| Missing ambiguity policy | 397.10 ms | 236.68 ms | 1.68× |
| Fractional ambiguity policy | 320.91 ms | 229.65 ms | 1.40× |
| Literal ambiguity policy | 320.69 ms | 230.74 ms | 1.39× |

### Geometric means

- all cases (119 cases), 2.53×
- algorithm only (47 cases), 8.79×
- all end-to-end and CLI cases (72 cases), 1.12×
- whole CLI (13 cases), 1.41×
- Python API (4 cases), 0.98×

### Remaining runtime

The stage profile uses three warm in-process repetitions. Large gappy spent
131 ms reading FASTA and 150 ms creating the alignment matrix, out of
418 ms total. Its trim step used 77 ms. Large KPI/KPIC and combined workloads
still pay for parsing and matrix construction, as well as the character
frequency pass. Large KPIC spent 176 ms trimming and 152 ms writing, out of
622 ms. Closely integrating FASTA ingestion and matrix construction would
address costs shared by these modes if it can retain format and header
semantics.

Medium smart-gap spent 130 ms writing, out of 300 ms. Large gappyout spent
304 ms writing, out of 688 ms. More of the gappyout gain must come from
output materialization or writing. Entropy and composition each spent
about 175 ms in trimming, in addition to roughly 260–280 ms in reading and
matrix construction. Sharing or narrowing the character counts is the next
candidate there.

The medium heterotachy guide-tree stage fell from 452 ms to 11 ms. Clade
scoring consumed 381 ms in the candidate, making it the next optimization
target. C3 spent about 3.2 of 3.8 ms reading and writing. Medium ECOMP spent
132 of 301 ms writing. Compression tuning must account for the encoded
length stored in the archive metadata. Cold startup still includes loading
NumPy and Biopython, and the API end-to-end sample remains close to parity.
Smaller inputs expose file access and CLI startup costs that the algorithm
benchmarks exclude.

The suite-wide end-to-end result remains below 1.75×. The largest
verified end-to-end gains are workload-specific, especially the mixed-case
paths. The tiny first-call heterotachy regression and the remaining
format-level timing differences are recorded in the interleaved audit.

The [54-case core checkpoint](benchmark-results/clipkit-175x-core-161eb0c.json)
also matched the reference. Its geometric means were 9.00× for algorithm-only
cases and 1.08× for end-to-end cases. The small core CLI subset averaged 1.21×.

### Heterotachy scoring follow-up

Commit `5a87a86` batches sites with the same observed-state count, preserving
the scalar reduction order and using the same gap and ambiguity filters.
The deterministic medium nucleotide CLI fixture produced the same FASTA
bytes across all three policies. Five interleaved measurements per policy
gave median wall times of 630 to 250 ms for missing, 690 to 250 ms for
fractional, and 450 to 240 ms for literal. Median peak RSS was 44,544 KiB on
both revisions. The samples are in the targeted benchmark JSON. Differential
tests match the rounded scalar scores on nucleotide and amino-acid inputs
that cross a column batch boundary. The full Linux test matrix passed
526 unit and 135 integration tests per interpreter at `5a87a86`.

### Incremental implementation history

Each row links the runtime change to its targeted verification. The targeted
JSON identifies the timing reference and candidate for each comparison.

| Commit | Change | Verification |
| --- | --- | --- |
| `eb97f94` | Defer Biopython tree construction imports | Exact help and version text with interleaved cold startup |
| `a63f9d1` | Infer alphabets from distinct symbols | Construction and gappy comparisons plus mixed-case and later-record tests |
| `84f582b` | Count small gap alphabets directly | Exact large gap scores and row-batch coverage |
| `45d439f` | Share the frequency table across report tracks | Exact HTML and FASTA with one frequency pass |
| `a710fc1` | Reuse normalized counts in mixed combined modes | Exact KPIC-gappy and KPIC smart-gap FASTA |
| `9eb38d6` | Normalize large ASCII arrays in bounded batches | Exact mixed-case trimming with Unicode fallback and strided-array tests |
| `e2eb3fd` | Benchmark startup, reports, and ambiguity handling | Exact comparisons for the added cases |
| `6b39356` | Reuse entropy totals for missing-site filtering | Exact entropy FASTA and analyzable-site fractions |
| `9b72df9` | Read standard large FASTA in blocks | Exact records, headers, fallback behavior, and trimmed output |
| `c3a70a0` | Batch Fitch scoring during the existing tree search | Exact topology, branch lengths, order, and heterotachy FASTA |
| `bb9de1a` | Retain raw timing samples and add larger workflows | Expected fixture shapes and exact added comparisons |
| `ab5fa77` | Use fast residue normalization in gap-only log checks | Exact debug log, FASTA, and first warned record |
| `de03988` | Audit CLI logs, complements, and tiny heterotachy | Exact combined artifacts and CLI output |
| `4b24f75` | Accumulate the ECOMP alphabet by sequence | Exact decoded alignment, metadata, and archive length |
| `684e8bb` | Batch wide FASTA line wrapping | Byte-identical output and complement, including header fallback |
| `5ed743a` | Fill very wide ASCII alignment matrices from bytes | Exact case, editable Unicode storage, and matching output |
| `8f3ecd2` | Suppress formatting work for disabled INFO messages | Exact quiet artifacts and entrypoint/state tests |
| `161eb0c` | Combine Boolean masks before converting to trim indices | Exact positions, thresholds, and classifications |
| `5a87a86` | Batch clade entropy scoring | Exact rounded scores and all ambiguity-policy FASTA outputs |

When reverting `45d439f`, revert `a710fc1` first because the mixed-case
schedule depends on normalized frequency reuse. Revert `ab5fa77` and
`5a87a86` before removing the fast normalization helper from `9eb38d6`.
The benchmark-only commits can be reverted independently of the production
code. The checkpoint and audit files identify their runtime revision.

### Wide FASTA and Python API follow-up

The targeted JSON also records the comparisons below. These are incremental
measurements against the runtime revision immediately before each change.

| Commit | Workflow | Reference median | Candidate median | Peak RSS |
| --- | --- | ---: | ---: | --- |
| `d9385d5` | Large CST with every site kept | 920 ms | 770 ms | 520,268 to 347,936 KiB |
| `d9385d5` | The same CST job with logging | 2.52 s | 2.17 s | 521,932 to 349,260 KiB |
| `a6edb03` | Small raw gappy API | 4.90 ms | 4.44 ms | Recorded in JSON |
| `a6edb03` | Small raw KPIC smart-gap API | 4.92 ms | 4.55 ms | Recorded in JSON |
| `a6edb03` | Small path-based gappy API | 4.44 ms | 4.02 ms | Recorded in JSON |
| `7d08ec7` | Large gappyout FASTA | 767 ms | 593 ms | 537 to 314 MB |
| `7d08ec7` | Large smart-gap FASTA | 783 ms | 607 ms | 536 to 314 MB |
| `7d08ec7` | Large KPIC FASTA | 690 ms | 617 ms | 375 to 312 MB |
| `7d08ec7` | Medium smart-gap FASTA | 353 ms | 295 ms | 261 to 167 MB |
| `7d08ec7` | Medium log and complement | 1.07 s | 1.01 s | 257 to 174 MB |

The keep-all path preserves an independent Bio alignment without copying
the input residue matrix first. Its two timing cases used the large
`EOG092C0CZK_aa_aln.fasta` fixture and a CST file that marks all 29,838 sites
`keep`. The API reads explicit ASCII FASTA from an in-memory text stream
with the same newline normalization as file input. The tests compare
headers, sequences, and statistics from CRLF text with the file path.

The wide FASTA exporter validates eligible matrices and positions before
opening the destination, then feeds one record at a time to the existing
wrapped writer. The legacy output conversion still handles other inputs
and formats. FASTA, log, and companion files matched exactly in the
streaming benchmark. A malformed-record regression test verifies that an
existing destination survives a conversion failure. Revert `7d08ec7`
before reverting the batched writer in `684e8bb`.

### Stop-codon masking follow-up

Commit `b78e9a0` checks gap membership on integer residue codes in bounded
row batches. The source matrix stays read-only until the selected stop
codons are masked. The targeted JSON contains five interleaved samples for
each case on the deterministic stop-codon fixture, with matching sequences,
mask counts, and emitted FASTA.

| Case | Reference median | Candidate median | Speedup | Peak RSS |
| --- | ---: | ---: | ---: | --- |
| Terminal masking, algorithm | 167 ms | 70 ms | 2.40× | 75.4 to 67.6 MB |
| Internal masking, end to end | 178 ms | 84 ms | 2.13× | 76.0 to 69.0 MB |
| All masking, end to end | 179 ms | 84 ms | 2.12× | 76.6 to 69.0 MB |

The unbounded prototype used about 20 MB more peak RSS than the reference.
The row batches retain the runtime gain and lower memory use. The focused
stop-codon, ambiguity-handling, and MSA suite passed 241 tests.

### CLI startup follow-up

Commit `f19e096` reads help and version arguments through the existing
parser before importing NumPy or Biopython. The console script and
`python -m clipkit` use that parser. Importing the callable from
`clipkit.api` still loads the trimming runtime, so the direct API retains
its existing first-call behavior. The file-format and sequence-type
enumerations retain their original import names for pickling.

| Case | Reference median | Candidate median | Speedup | Peak RSS |
| --- | ---: | ---: | ---: | --- |
| CLI help | 215 ms | 46 ms | 4.69× | 44.0 to 22.0 MB |
| CLI version | 224 ms | 47 ms | 4.81× | 44.0 to 22.0 MB |
| Small gappy CLI | 230 ms | 225 ms | 1.02× | 47.2 to 46.4 MB |
| Medium nucleotide smart-gap CLI | 245 ms | 237 ms | 1.03× | 48.8 MB on both |
| Small raw gappy API | 4.34 ms | 4.32 ms | 1.01× | 48.8 MB on both |

The interleaved comparison kept six fresh-process samples per source and
checked exact outputs. The 25 focused tests verify the public callable,
enum pickle round trips, and startup imports. The installed script's help
and version text also matched the reference module invocation.

### Matrix construction and counting follow-up

Commit `eabf605` packs aligned ASCII Bio.Seq bytes before widening them in
one NumPy operation. The returned matrix retains its writable Unicode
representation, case, and inferred sequence type. Tests also check the
original error path for non-ASCII residues. Commit `9908ef8` uses smaller
batches for character histograms, based on measurements on sparse
amino-acid and dense amino-acid and nucleotide alignments. Heterotachy
retains its independent site batch width.

| Change | Workload | Reference median | Candidate median |
| --- | --- | ---: | ---: |
| Packed input | Large amino-acid matrix construction | 159 ms | 66 ms |
| Packed input | Medium amino-acid matrix construction | 80 ms | 33 ms |
| Packed input | Large gappy end to end | 476 ms | 349 ms |
| Packed input | Medium smart-gap end to end | 292 ms | 232 ms |
| Smaller count batches | Large KPIC end to end | 504 ms | 474 ms |
| Smaller count batches | Large entropy end to end | 520 ms | 493 ms |
| Smaller count batches | Large composition bias end to end | 537 ms | 505 ms |
| Smaller count batches | Dense amino-acid KPI algorithm | 20.1 ms | 13.4 ms |

Five interleaved samples per source matched the case-specific output and
metadata. Peak RSS for the listed large and medium workflows stayed within
1%, or fell. The [instrumented stage samples](benchmark-results/clipkit-175x-stage-profile.json)
show large gappy parsing and matrix construction at 134 ms and 67 ms,
respectively, on `9908ef8`. Large KPIC spent 128 ms trimming and 89 ms
writing. The medium smart-gap writer consumed 84 of 213 ms, and the
medium ECOMP writer consumed 134 of 261 ms. Medium heterotachy clade
scoring fell to 3 ms after batching.

The exact 119-case comparison for `f19e096` measured a 1.29× geometric
mean across end-to-end workloads and 2.00× for the whole-CLI subset.
That comparison precedes the packed-input and histogram changes.

### Direct FASTA row output

Commit `ae894a4` converts selected ASCII code points directly to bytes
and shares the line-wrap writer with the Bio path. Streaming now covers
outputs starting at 2,048 columns. The destination remains untouched
when unequal decoded row lengths would fail Bio validation, including
matrices with trailing NUL padding.

| Case | Reference median | Candidate median |
| --- | ---: | ---: |
| Medium C3 nucleotide | 6.38 ms | 5.31 ms |
| Large gappyout | 485 ms | 466 ms |
| Large KPIC | 483 ms | 465 ms |
| Medium smart-gap | 231 ms | 223 ms |
| Medium log and complement | 939 ms | 925 ms |

Six interleaved samples per tree matched exact output and companion
files. Peak RSS was flat within about 1%. The input-parser threshold
stayed at 1 MiB after smaller blocks and earlier eligibility had mixed
results across medium input shapes.

### Measured checkpoint at `9908ef8`

The [119-case comprehensive results](benchmark-results/clipkit-175x-comprehensive-9908ef8.json)
compare `9908ef8` with `6044b1e` using the recorded CPython 3.12.14,
NumPy 2.0.2, and Biopython 1.83 environment on the same Linux host. They
contain five fresh-process measurements per case and source and record
exact output agreement. This checkpoint includes packed input and
smaller histograms. The direct FASTA row follow-up has its own targeted
comparison above.

### Same-machine suite comparison

The comprehensive suite compared 119 cases exactly against `6044b1e`.
Each row shows a representative workload for one public trimming mode.

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| `gappy` | 623.25 ms | 356.29 ms | 1.75× |
| `smart_gap` | 1,013.93 ms | 484.38 ms | 2.09× |
| `block_gappy` | 646.82 ms | 364.80 ms | 1.77× |
| `gappyout` | 1,023.80 ms | 484.44 ms | 2.11× |
| `kpi` | 755.24 ms | 463.43 ms | 1.63× |
| `kpic` | 784.87 ms | 468.31 ms | 1.68× |
| `kpi_gappy` | 646.98 ms | 424.00 ms | 1.53× |
| `kpic_gappy` | 631.18 ms | 423.81 ms | 1.49× |
| `kpi_smart_gap` | 738.11 ms | 470.06 ms | 1.57× |
| `kpic_smart_gap` | 756.53 ms | 484.36 ms | 1.56× |
| `entropy` | 1,005.06 ms | 497.72 ms | 2.02× |
| `composition_bias` | 1,035.69 ms | 505.03 ms | 2.05× |
| `heterotachy` | 1,084.93 ms | 247.77 ms | 4.38× |
| `cst` | 0.69 ms | 0.59 ms | 1.18× |
| `c3` | 6.35 ms | 6.20 ms | 1.02× |

### Major workflows

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| CLI help startup | 237.00 ms | 45.41 ms | 5.22× |
| CLI version startup | 233.53 ms | 45.88 ms | 5.09× |
| CLI small FASTA trimming | 251.82 ms | 220.39 ms | 1.14× |
| Python API from file | 4.38 ms | 3.93 ms | 1.12× |
| Python API from string | 4.71 ms | 4.20 ms | 1.12× |
| HTML report | 1,820.13 ms | 778.74 ms | 2.34× |
| Log and complementary output | 1,190.12 ms | 941.97 ms | 1.26× |
| Mixed-case amino-acid CLI | 1,050.96 ms | 253.60 ms | 4.14× |
| ECOMP input | 2.19 ms | 2.00 ms | 1.09× |
| ECOMP output | 403.95 ms | 282.75 ms | 1.43× |
| Codon trimming | 5.22 ms | 5.00 ms | 1.05× |
| Terminal stop-codon masking | 185.19 ms | 75.46 ms | 2.45× |
| Internal stop-codon masking | 181.15 ms | 77.81 ms | 2.33× |
| All stop-codon masking | 180.90 ms | 77.25 ms | 2.34× |
| One requested thread | 463.95 ms | 231.33 ms | 2.01× |
| Four requested threads | 445.84 ms | 234.18 ms | 1.90× |
| Missing ambiguity policy | 397.10 ms | 218.78 ms | 1.82× |
| Fractional ambiguity policy | 320.91 ms | 217.70 ms | 1.47× |
| Literal ambiguity policy | 320.69 ms | 218.04 ms | 1.47× |

### Geometric means

- all cases (119 cases), 3.43×
- algorithm only (47 cases), 12.96×
- all end-to-end and CLI cases (72 cases), 1.44×
- whole CLI (13 cases), 2.04×
- Python API (4 cases), 1.23×

The representative large gap modes, entropy, composition bias, and
medium heterotachy reached approximately 1.75× or better at this checkpoint.
The combined gap and classification modes reached 1.49–1.57×, and the
tiny CST and medium C3 cases reached 1.18× and 1.02×. The aggregate
end-to-end result of 1.44× leaves the suite-wide target open. The stage
samples show the largest shared costs still include FASTA reading,
residue counting, and output serialization.

The Linux test matrix for the later runtime revision `ae894a4` passed
533 unit and 135 integration tests on each of CPython 3.12.14, 3.13.15,
3.14.7, and 3.15.2.

### Sparse and medium residue counts

Commit `93a6e74` applies a single residue filter to configured gaps and
excluded ambiguity symbols. Its interleaved timing and memory medians
stayed within 1%. Commit `218368e` bins the non-dash residues in highly
sparse batches and reconstructs exact dash counts by subtraction. A
sample selects likely sparse input, while a per-batch density check
retains the full histogram for dense regions. Both dense controls
stayed within 2% in time and 3% in peak RSS.

Commit `c00fc29` reuses the bounded integer membership pass for
medium-sized ASCII gap alphabets, including mixed-case nucleotide
input after normalization. Tests compare the rounded result with
Unicode membership, including an uppercase expansion.

| Change | Workload | Reference median | Candidate median |
| --- | --- | ---: | ---: |
| Sparse histogram | Large KPIC algorithm | 207 ms | 187 ms |
| Sparse histogram | Large KPIC end to end | 471 ms | 451 ms |
| Sparse histogram | Large KPIC-gappy end to end | 421 ms | 399 ms |
| Integer membership | Medium nucleotide ends-only | 16.14 ms | 8.27 ms |
| Integer membership | Medium nucleotide gappy algorithm | 15.31 ms | 7.32 ms |
| Integer membership | Small raw gappy API | 4.25 ms | 2.94 ms |

The samples in the targeted JSON match the checked outputs, positions,
and metadata. The preceding `ae894a4` checkpoint compared all 119
comprehensive cases exactly and reached 1.45× across end-to-end
workloads, 2.04× for whole CLI runs, and 1.21× for Python API cases.
The focused changes above follow that checkpoint.

### Stop-codon selection and cold execution audit

Commit `4f44170` finds the last complete codon in each row with one
indexed operation and reuses the stop mask for internal and all-mode
selection. The focused unit tests include a row composed entirely of
gaps. Six fresh-process observations per side were interleaved after
two warm-ups for all three stop modes.

| Stop-codon mode | Previous end-to-end median | Indexed end-to-end median | Speedup |
| --- | ---: | ---: | ---: |
| Terminal | 74.56 ms | 71.58 ms | 1.04× |
| Internal | 75.20 ms | 71.62 ms | 1.05× |
| All | 75.64 ms | 72.77 ms | 1.04× |

Outputs and masking metadata matched exactly. Median peak RSS fell
1.4–2.0 MiB for terminal and internal modes and rose 0.3 MiB for all-mode.

The interleaved audit of `c00fc29` against `6044b1e` compared 18 cases
exactly. Medium nucleotide ends-only went from 89.67 to 7.89 ms, and large
KPIC-smart-gap went from 754.87 to 442.72 ms. Tiny heterotachy went from
6.06 to 25.55 ms in the first execution after importing ClipKIT because
the tree-construction import now falls inside that timer. A fresh CLI
invocation of the same workflow went from 246.22 to 239.71 ms.

### Updated 119-case checkpoint

The [raw comparison](benchmark-results/clipkit-175x-comprehensive-c00fc29.json)
records the five fresh-process measurements per case and revision at
`c00fc29`. The representative combined modes improved to 1.59–1.68×.

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| `gappy` | 623.25 ms | 361.17 ms | 1.73× |
| `smart_gap` | 1,013.93 ms | 466.87 ms | 2.17× |
| `block_gappy` | 646.82 ms | 371.49 ms | 1.74× |
| `gappyout` | 1,023.80 ms | 464.61 ms | 2.20× |
| `kpi` | 755.24 ms | 445.83 ms | 1.69× |
| `kpic` | 784.87 ms | 454.53 ms | 1.73× |
| `kpi_gappy` | 646.98 ms | 397.75 ms | 1.63× |
| `kpic_gappy` | 631.18 ms | 397.65 ms | 1.59× |
| `kpi_smart_gap` | 738.11 ms | 442.95 ms | 1.67× |
| `kpic_smart_gap` | 756.53 ms | 450.47 ms | 1.68× |
| `entropy` | 1,005.06 ms | 468.71 ms | 2.14× |
| `composition_bias` | 1,035.69 ms | 472.42 ms | 2.19× |
| `heterotachy` | 1,084.93 ms | 246.13 ms | 4.41× |
| `cst` | 0.69 ms | 0.70 ms | 0.99× |
| `c3` | 6.35 ms | 5.38 ms | 1.18× |

### Major workflows

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| CLI help startup | 237.00 ms | 46.92 ms | 5.05× |
| CLI version startup | 233.53 ms | 45.77 ms | 5.10× |
| CLI small FASTA trimming | 251.82 ms | 223.74 ms | 1.13× |
| Python API from file | 4.38 ms | 2.60 ms | 1.69× |
| Python API from string | 4.71 ms | 2.87 ms | 1.64× |
| HTML report | 1,820.13 ms | 746.51 ms | 2.44× |
| Log and complementary output | 1,190.12 ms | 911.32 ms | 1.31× |
| Mixed-case amino-acid CLI | 1,050.96 ms | 248.10 ms | 4.24× |
| ECOMP input | 2.19 ms | 1.58 ms | 1.38× |
| ECOMP output | 403.95 ms | 280.31 ms | 1.44× |
| Codon trimming | 5.22 ms | 2.13 ms | 2.46× |
| Terminal stop-codon masking | 185.19 ms | 76.61 ms | 2.42× |
| Internal stop-codon masking | 181.15 ms | 75.65 ms | 2.39× |
| All stop-codon masking | 180.90 ms | 75.79 ms | 2.39× |
| One requested thread | 463.95 ms | 225.98 ms | 2.05× |
| Four requested threads | 445.84 ms | 225.24 ms | 1.98× |
| Missing ambiguity policy | 397.10 ms | 220.15 ms | 1.80× |
| Fractional ambiguity policy | 320.91 ms | 218.99 ms | 1.47× |
| Literal ambiguity policy | 320.69 ms | 220.23 ms | 1.46× |

### Geometric means

- all cases (119 cases), 3.66×
- algorithm only (47 cases), 13.88×
- all end-to-end and CLI cases (72 cases), 1.53×
- whole CLI (13 cases), 2.05×
- Python API (4 cases), 1.48×

### Sequence-type inference with a gappy first row

Commit `1ddf348` uses a literal deletion table for single-character gaps.
Inference still checks additional records when the first row contains
fewer than 200 non-gap residues. On the medium amino-acid fixture with
its first row replaced by dashes, six interleaved observations per side
measured 1,750.28 to 298.15 ms (5.87×) and median peak RSS fell from 233
to 159 MiB. The unmodified medium fixture measured 222.01 to 221.16 ms,
and the small gappy fixture measured 4.34 to 4.12 ms. All three cases
matched their reference outputs and metadata exactly.

The [updated 54-case core checkpoint](benchmark-results/clipkit-175x-core-c00fc29.json)
also matched `6044b1e` exactly. Its geometric means were 1.40× for
end-to-end and CLI cases, 1.25× for whole CLI invocation, 1.29× for the
Python API subset, and 13.16× for algorithm-only cases. The
comprehensive checkpoint included a 6.9% RSS increase for the small
raw Python API case. The interleaved audit measured the same median
RSS for both revisions on that case.

### Gap-only warnings for log output

Commit `4ed4dc8` tests wide ASCII rows for resolved residues with byte
deletion, including matching lowercase gap symbols. Other rows retain
Unicode normalization and membership. The warning still reports the
first matching trimmed row. A focused test covers Unicode characters
that map into the ASCII gap alphabet when uppercased.

The warning stage fell from about 173 to 10 ms on the medium log-and-
complement sample. Six interleaved observations per revision measured
1.22× for the whole command, with exact output and a 0.24% median peak
RSS increase. The small control varied by 2.7% in time. The remaining
cost in the medium workflow is debug log generation and emission,
which consumed about 279 ms after the change.

Commit `c9774e8` shares the classification count table with gappyness
when writing a per-site log. Six interleaved measurements per revision
showed another 1.05× improvement for the whole command, with exact
output and a 0.24% increase in median peak RSS. An interleaved
comparison with the original baseline measured 1,230.88 to 722.26 ms
(1.70×) for log and complement, and 243.12 to 218.52 ms (1.11×) for
the small standalone CLI control.

On Linux x86-64, `make test` on `c9774e8` passed 538 unit and 135
integration tests on each of CPython 3.12.14, 3.13.15, 3.14.7, and
3.15.2. Each interpreter used its own temporary directory.

### Checkpoint after sequence-type inference optimization

The [core comparison](benchmark-results/clipkit-175x-core-0b9ceed.json) and
[comprehensive comparison](benchmark-results/clipkit-175x-comprehensive-0b9ceed.json)
for `0b9ceed` matched all 54 and 119 cases, respectively, against `6044b1e`.

| Geometric mean | Core | Comprehensive |
| --- | ---: | ---: |
| End-to-end and CLI | 1.42× | 1.54× |
| Whole CLI | 1.31× | 2.03× |
| Python API | 1.35× | 1.48× |
| Algorithm only | 13.33× | 13.95× |

The comprehensive run flagged tiny heterotachy, FASTA output, and Mauve
output above 5% slower. The earlier interleaved audit reproduced the
first-call heterotachy deferral, while FASTA and Mauve were faster in
that audit. The log-path improvements follow this checkpoint.

### Gap membership on highly sparse alignments

Commit `8fbeef9` reuses the existing dash-density sample and counts
resolved residues in sparse row batches. Subtracting those counts from
the batch size gives exact gap and ambiguity counts without traversing
the full matrix for each gap symbol. Dense batches use the existing
membership pass. The new regression test combines sparse rows,
non-dash ambiguity symbols, and a dense batch in one alignment.

| Workload | Previous median | Candidate median | Speedup |
| --- | ---: | ---: | ---: |
| Medium smart-gap | 219.26 ms | 207.38 ms | 1.06× |
| Large gappy end to end | 354.41 ms | 344.49 ms | 1.03× |
| Large gappy algorithm | 140.45 ms | 139.33 ms | 1.01× |
| Dense amino-acid gappy algorithm | 10.04 ms | 10.11 ms | 0.99× |

The six interleaved observations per side matched exactly. Median peak
RSS stayed within 1% across these workloads.

An interleaved comparison with `6044b1e` measured 632.53 to 344.83 ms
(1.83×) for large gappy, and 464.39 to 208.32 ms (2.23×) for
medium smart-gap. Both matched the reference exactly. The Linux test
matrix passed 539 unit and 135 integration tests on each of CPython
3.12.14, 3.13.15, 3.14.7, and 3.15.2 at `8fbeef9`.

### Keep-only custom site trimming

Commit `7821042` builds a boolean mask from validated keep indices and
extracts the complementary positions in order. Commit `c262934` passes
that array through trimming and validation. The public list helper
still returns the same values and element types. Revert `c262934`
before reverting `7821042` to restore the former list-based execution.

The benchmark alignment contains four identical rows made by repeating
`ACGT` 250,000 times. The CST file retains one-based positions 1, 38,
75, and every 37th position thereafter. Each comparison interleaved six
fresh-process samples per side following two warm-ups.

| Revision | End-to-end median | Median peak RSS |
| --- | ---: | ---: |
| Original baseline `6044b1e` | 373.82 ms | 126.9 MiB |
| Keep-only array `c262934` | 68.88 ms | 88.7 MiB |

Trimmed FASTA matched exactly. The 5.43× gain applies to this wide
keep-only workload. The tiny CST case with an explicit trim directive
measured 0.70 to 0.63 ms in the baseline comparison.

The [119-case log-optimization checkpoint](benchmark-results/clipkit-175x-comprehensive-c9774e8.json)
matched `6044b1e` exactly before the sparse gap and CST follow-ups. It
measured 1.56× across end-to-end cases, 2.10× for the whole-CLI subset,
1.50× for Python API cases, and 13.91× for algorithm-only cases. Tiny
first-call heterotachy was the only runtime regression above 5% in
that suite.

### C3 and small-file profiling

At `c262934`, median warmed C3 execution took 3.53 ms. Reading took
1.23 ms and writing took 1.76 ms, together accounting for 85% of the
profiled time. Trimming itself took 0.05 ms. Medium ECOMP spent 133 ms
of 244 ms in writing, which includes compression. The stage JSON
records both profiles and the small gappy control.

A prototype gathered all selected columns into one array for medium
FASTA output. It slowed fresh-process C3 by 8.7% and codon trimming by
7.8%, so the implementation retains row-wise selection. Commit
`9bce4f6` removes duplicate Biopython alignment-length traversals.
The interleaved C3 median moved from 5.35 to 5.25 ms with exact output.
Further gains on this C3 case depend on parsing and writing costs.

The [54-case core result](benchmark-results/clipkit-175x-core-c9774e8.json)
at `c9774e8` matched the baseline exactly. The end-to-end geometric mean
was 1.42× and the algorithm-only geometric mean was 13.25×.

### Remaining time to reach 1.75×

The values below use the same-machine `c9774e8` comprehensive checkpoint.
The last column is the additional reduction required from the measured
candidate median to reach 1.75× for that case.

| Workflow | Baseline | Candidate | Speedup | Additional reduction |
| --- | ---: | ---: | ---: | ---: |
| Medium C3 | 6.35 ms | 5.24 ms | 1.21× | 1.62 ms |
| Tiny explicit CST | 0.69 ms | 0.63 ms | 1.10× | 0.23 ms |
| Large KPI-gappy | 646.98 ms | 390.69 ms | 1.66× | 20.99 ms |
| Large KPIC-gappy | 631.18 ms | 398.91 ms | 1.58× | 38.23 ms |
| Large KPI-smart-gap | 738.11 ms | 423.03 ms | 1.74× | 1.25 ms |
| Large KPIC-smart-gap | 756.53 ms | 449.46 ms | 1.68× | 17.16 ms |
| Small gappy CLI | 251.82 ms | 211.43 ms | 1.19× | 67.53 ms |
| Small raw gappy API | 4.71 ms | 2.87 ms | 1.64× | 0.18 ms |
| ECOMP input | 2.19 ms | 1.40 ms | 1.56× | 0.15 ms |
| Medium ECOMP output | 403.95 ms | 279.56 ms | 1.44× | 48.73 ms |
| Medium fractional-ambiguity CLI | 320.91 ms | 216.41 ms | 1.48× | 33.03 ms |
| Medium literal-ambiguity CLI | 320.69 ms | 217.39 ms | 1.48× | 34.14 ms |
| Medium log and complement CLI | 1,190.12 ms | 709.96 ms | 1.68× | 29.89 ms |

C3 needs faster FASTA parsing or serialization. Tiny CST also has little
algorithm time to remove, while the separate million-site CST workload
benefited from retaining site arrays. The large combined modes share the
residue histogram already. Further gains must come from faster histogram
construction, parsing, or encoding the selected rows. Medium ECOMP output
needs a reduction in serialization or gzip time while preserving decoded
content and archive metadata. The small CLI and both ambiguity CLI cases
spend substantially more time starting the process than trimming, so measure
imports and the full command together when evaluating deferred imports.
The small in-process API and ECOMP reader have less than 0.2 ms each to
recover. The log workflow still spends most of its post-trim time building
and emitting the per-site log.

Two follow-up prototypes were evaluated with six interleaved
fresh-process observations per source. Taking the alignment length from
the first parsed row left the four controls within 2%, so the shared
Biopython length call remains. Deferring `Bio.SeqIO` until parsing or
writing moved the small in-process FASTA call from 4.14 to 24.15 ms
while changing the full small and medium log CLI invocations by less
than 1%. The existing import remains. Both prototypes matched the
checked outputs exactly.

### ECOMP fallback input and output

Commit `3a87d78` decodes decompressed FASTA payloads of at least 1 MiB
incrementally. The parser continues to split on LF and strip line-edge
whitespace. Tests cover Unicode headers, CRLF, and a final line without
a terminator on both size paths.

The large ECOMP fixture wraps the amino-acid alignment
`tests/integration/samples/EOG092C0CZK_aa_aln.fasta` in the existing gzip
fallback archive format. Six interleaved observations per side measured
559.34 to 358.34 ms (1.56×) compared to `9bce4f6`, and 815.33 to
361.17 ms (2.26×) compared to `6044b1e`. The incremental comparison
reduced median peak RSS from 303.2 to 300.1 MiB. The small control
stayed within 2%, and the checked outputs matched exactly.

The medium ECOMP output profile measured about 113 ms in gzip
compression, 7 ms preparing metadata, and 3 ms encoding metadata.
Sharing the converted sequence strings moved the whole median from
263.14 to 262.74 ms. The existing output path remains.

The [119-case checkpoint at `c262934`](benchmark-results/clipkit-175x-comprehensive-c262934.json)
matched `6044b1e` exactly. The geometric means were 1.58× end to end,
2.12× for the whole CLI subset, 1.53× for Python API cases, and 13.95×
for algorithm-only cases. It precedes incremental ECOMP decoding.

### Buffered per-site log file writes

Commit `38d8946` buffers writes to ClipKIT's per-site log file and
flushes after the log records are written. Other handlers continue to
receive each formatted logging event through the logging framework.
The file handler is closed on exit and logger settings are restored
even if closing raises.

Six interleaved observations reduced the medium log-and-complement
median from 718.24 to 639.92 ms against `3a87d78`, with peak RSS
essentially equal. A separate interleaved comparison against `6044b1e`
measured 1,242.04 to 663.39 ms (1.87×). The log, complement, and main
alignment matched exactly. The small no-log control stayed within 1%
in the incremental comparison.

On `38d8946`, Linux `make test` passed 540 unit and 135 integration
tests on each of CPython 3.12.14, 3.13.15, 3.14.7, and 3.15.2.

The earlier [core comparison at `c262934`](benchmark-results/clipkit-175x-core-c262934.json)
matched all 54 cases exactly and measured 1.42× end to end and 13.26×
for algorithm-only cases.

### Direct row-byte matrix construction

Commit `a9ccb80` widens ASCII bytes from each row directly into the
final writable Unicode matrix. This eliminates the temporary packed
copy of the entire alignment. A focused test checks that shorter and
longer records modified after Biopython validation still produce the
previous padding and truncation, including embedded NUL bytes.

| Workload | Previous median | Candidate median | Peak RSS before and after |
| --- | ---: | ---: | ---: |
| Large amino-acid construction | 66.69 ms | 63.24 ms | 300.1 → 258.1 MiB |
| Medium amino-acid construction | 32.27 ms | 31.29 ms | 156.6 → 138.3 MiB |
| Large KPIC-gappy end to end | 402.05 ms | 387.01 ms | 302.9 → 277.8 MiB |
| Dense nucleotide KPIC algorithm | 62.42 ms | 63.04 ms | 73.2 → 70.2 MiB |

Each comparison used six interleaved fresh-process observations per
source. Output and the case-specific correctness metadata matched
exactly.

### Single-byte FASTA header scan

Commit `192f57a` splits standard FASTA blocks on the one-byte header
marker and retains newline-aware splitting for markers within a title
or sequence. Tests compare the records for both cases with Biopython.
Six interleaved measurements per side improved large gappy from 341.73
to 288.80 ms, medium smart-gap from 205.60 to 181.11 ms, and large
KPIC-gappy from 393.53 to 337.01 ms. Peak RSS stayed within 1%.

| Large combined mode | Baseline median | Candidate median | Speedup |
| --- | ---: | ---: | ---: |
| KPI-gappy | 659.12 ms | 377.50 ms | 1.75× |
| KPIC-gappy | 628.46 ms | 347.65 ms | 1.81× |
| KPI-smart-gap | 728.19 ms | 376.44 ms | 1.93× |
| KPIC-smart-gap | 759.45 ms | 386.95 ms | 1.96× |

The baseline comparison interleaved six fresh-process observations per
side and matched the checked outputs exactly. Commit `913afff` then
extended the same reader to regular files starting at 64 KiB. Small
amino-acid trimming improved from 4.14 to 3.81 ms and the file API
from 2.50 to 2.28 ms, with exact outputs. C3 and codon controls stayed
within 3%. The smaller cutoff depends on the preceding header scan
improvement for its timing.

The [comprehensive checkpoint at `3a87d78`](benchmark-results/clipkit-175x-comprehensive-3a87d78.json)
matched all 119 cases against `6044b1e`. It measured 1.59× end to end,
2.11× for whole CLI runs, 1.50× for Python API cases, and 14.02× for
algorithm-only cases. It precedes the FASTA header scan changes.

### Block parsing for the raw ASCII FASTA API

Commit `2ab1d58` sends raw ASCII FASTA strings of at least 64 KiB through
the existing block reader using an in-memory byte stream. Legacy text
handling remains available when the block reader declines the input.
Tests compare trimmed records and headers for LF and CRLF inputs and check
that neither path writes a temporary alignment file. This optimization
uses the block reader and cutoff developed in `192f57a` and `913afff`.

| Small public API workflow | Baseline median | Candidate median | Speedup |
| --- | ---: | ---: | ---: |
| Raw gappy | 4.84 ms | 2.14 ms | 2.26× |
| Raw KPIC-smart-gap | 5.37 ms | 3.51 ms | 1.53× |
| File gappy | 4.44 ms | 2.28 ms | 1.95× |

Six interleaved observations per source agreed exactly on the checked
artifacts and metadata. Compared with the preceding runtime revision,
raw gappy fell from 2.77 to 2.09 ms and raw KPIC-smart-gap fell from 3.92
to 3.34 ms. Median peak RSS was flat or lower in those two comparisons.

On Linux, `make test` at `2ab1d58` passed 545 unit and 135 integration
tests on each of CPython 3.12.14, 3.13.15, 3.14.7, and 3.15.2. The earlier
[core checkpoint at `3a87d78`](benchmark-results/clipkit-175x-core-3a87d78.json)
matched all 54 cases and measured 1.43× end to end and 13.32× for
algorithm-only cases.

### Smart-gap slopes and small residue-count batches

Commit `a9e75d1` computes the existing slope formula in NumPy.
It retains the cumulative sum, transition order, duplicate filtering,
and first-half selection. Large KPIC-smart-gap improved from 391.62 to
385.05 ms. In-process timing for 1,000 points fell from about 691 to 46 μs.

Commit `fc041aa` halves the histogram column batch for small matrices
with compact residue alphabets. The small raw KPIC-smart-gap API
improved from 3.33 to 3.07 ms in the incremental comparison. A wider
batch slowed that case and the medium nucleotide algorithm case. The
interleaved baseline comparison measured 5.06 to 3.19 ms (1.59×) for
raw KPIC-smart-gap and 5.29 to 2.19 ms (2.42×) for raw gappy. Outputs
matched exactly. The [raw samples](benchmark-results/clipkit-175x-targeted.json)
include the incremental, baseline, and discarded prototype timings.

A separate interleaved audit measured 1.04× for tiny FASTA output and
1.07× for tiny Mauve output. It also measured 1.06× for a fresh tiny
heterotachy CLI process. The first-call in-process tree import remains
a separate cost. On Linux, `make test` at `fc041aa` passed 546 unit and
135 integration tests per interpreter on CPython 3.12.14, 3.13.15,
3.14.7, and 3.15.2.

### Comprehensive checkpoint after FASTA parsing changes

The [119-case comparison at `913afff`](benchmark-results/clipkit-175x-comprehensive-913afff.json)
matched `6044b1e` exactly. It precedes the raw API and smart-gap
follow-ups. Thirteen representative trimming-mode cases reached 1.75×
or better. The explicit tiny CST and medium C3 cases measured 1.10×
and 1.19×, respectively.

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| `gappy` | 623.25 ms | 290.67 ms | 2.14× |
| `smart_gap` | 1,013.93 ms | 390.03 ms | 2.60× |
| `block_gappy` | 646.82 ms | 283.62 ms | 2.28× |
| `gappyout` | 1,023.80 ms | 394.65 ms | 2.59× |
| `kpi` | 755.24 ms | 380.71 ms | 1.98× |
| `kpic` | 784.87 ms | 383.77 ms | 2.05× |
| `kpi_gappy` | 646.98 ms | 346.51 ms | 1.87× |
| `kpic_gappy` | 631.18 ms | 349.07 ms | 1.81× |
| `kpi_smart_gap` | 738.11 ms | 392.45 ms | 1.88× |
| `kpic_smart_gap` | 756.53 ms | 392.52 ms | 1.93× |
| `entropy` | 1,005.06 ms | 410.68 ms | 2.45× |
| `composition_bias` | 1,035.69 ms | 417.93 ms | 2.48× |
| `heterotachy` | 1,084.93 ms | 248.79 ms | 4.36× |
| `cst` | 0.69 ms | 0.63 ms | 1.10× |
| `c3` | 6.35 ms | 5.35 ms | 1.19× |

### Major workflows

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| CLI help startup | 237.00 ms | 46.68 ms | 5.08× |
| CLI version startup | 233.53 ms | 46.25 ms | 5.05× |
| CLI small FASTA trimming | 251.82 ms | 223.08 ms | 1.13× |
| Python API from file | 4.38 ms | 2.35 ms | 1.87× |
| Python API from string | 4.71 ms | 2.80 ms | 1.68× |
| HTML report | 1,820.13 ms | 674.87 ms | 2.70× |
| Log and complementary output | 1,190.12 ms | 618.36 ms | 1.92× |
| Mixed-case amino-acid CLI | 1,050.96 ms | 257.80 ms | 4.08× |
| ECOMP input | 2.19 ms | 1.40 ms | 1.56× |
| ECOMP output | 403.95 ms | 249.40 ms | 1.62× |
| Codon trimming | 5.22 ms | 2.01 ms | 2.60× |
| Terminal stop-codon masking | 185.19 ms | 85.50 ms | 2.17× |
| Internal stop-codon masking | 181.15 ms | 86.32 ms | 2.10× |
| All stop-codon masking | 180.90 ms | 85.83 ms | 2.11× |
| One requested thread | 463.95 ms | 188.11 ms | 2.47× |
| Four requested threads | 445.84 ms | 185.61 ms | 2.40× |
| Missing ambiguity policy | 397.10 ms | 227.64 ms | 1.74× |
| Fractional ambiguity policy | 320.91 ms | 231.79 ms | 1.38× |
| Literal ambiguity policy | 320.69 ms | 228.55 ms | 1.40× |

### Geometric means

- all cases (119 cases), 3.79×
- algorithm only (47 cases), 13.73×
- all end-to-end and CLI cases (72 cases), 1.64×
- whole CLI (13 cases), 2.10×
- Python API (4 cases), 1.55×

### Small combined-mode profiling and C3 fixture comparison

The warmed small raw KPIC-smart-gap API case at `fc041aa` spent about
0.64 ms parsing, 0.24 ms constructing the alignment matrix, and 0.50 ms
counting residues, out of 1.83 ms measured in process. Smart-gap
thresholding took about 0.11 ms, including 0.025 ms for slopes. A
column-major histogram layout left the measured controls within 3%
and did not improve either small combined case.

| C3 workload | Baseline median | Candidate median | Speedup |
| --- | ---: | ---: | ---: |
| Medium nucleotide | 6.55 ms | 5.48 ms | 1.19× |
| Codon fixture | 1.79 ms | 1.57 ms | 1.14× |
| Generated dense nucleotide (512 × 9,000) | 71.15 ms | 42.53 ms | 1.67× |

All three C3 cases matched the baseline exactly across six interleaved
observations per source. The larger case needs about another 1.9 ms
reduction to meet 1.75×. The earlier medium profile puts most time in
parsing and serialization.

The [54-case core comparison at `913afff`](benchmark-results/clipkit-175x-core-913afff.json)
also matched the baseline exactly.

### Direct FASTA output on smaller selections

Commit `fa671f8` applies the validated FASTA writer at 512 selected
columns. Small gappy fell from 3.98 to 3.69 ms in six interleaved
observations per source, with exact output and equal median peak RSS.
Keep and complement tests compare both selection sizes with Biopython.
A 128-column cutoff and a wrapped-indexing prototype offered too
little benefit to retain.

In a separate interleaved comparison with `6044b1e`, small gappy
measured 1.55×, small KPIC-gappy 1.32×, and medium log and complement
1.96×, with exact checked artifacts.

Commit `261b4b3` passes the existing wrapping array to the ASCII
decoder without materializing a separate bytes object. The median
for large KPIC-smart-gap fell from 386.48 to 375.14 ms against the
previous revision. Small gappy, medium C3, and medium smart-gap
stayed within 2%. All checked bytes matched.

### Cache-sized FASTA output validation

Commit `fef1ffa` validates each group of rows with adjacent minimum and
maximum reductions before opening the destination. The focused
validation test covers a non-ASCII residue in a later batch and checks
that an existing destination is preserved. Six interleaved observations
reduced large gappy from 286.19 to 279.04 ms and large KPIC-smart-gap
from 388.43 to 377.50 ms. Peak RSS and both small controls stayed
within 1%. Checked outputs matched exactly.

The earlier [119-case comparison at `fc041aa`](benchmark-results/clipkit-175x-comprehensive-fc041aa.json)
matched `6044b1e` exactly. Its geometric means were 1.66× end to end,
2.12× for the whole CLI subset, 1.76× for Python API cases, and 13.77×
for algorithm-only cases. Tiny first-call heterotachy was the only
runtime regression above 5%. The raw small gappy API case showed a
6.9% RSS increase in the comprehensive run and equal median RSS in
the interleaved audit.

### Reused FASTA line buffer

Commit `703e30d` allocates the direct FASTA wrapping buffer once for
each written selection. Its focused test compares two different
standard rows and an unusual-header fallback against Biopython.
Six interleaved observations per source measured 3.81 to 3.67 ms for
small gappy and 5.28 to 5.15 ms for medium C3. The large controls
stayed within 2%, and peak RSS was flat or lower except for a 0.4%
increase in large KPIC-smart-gap. Checked output bytes matched.

The [54-case core comparison at `fc041aa`](benchmark-results/clipkit-175x-core-fc041aa.json)
matched `6044b1e` exactly. It measured 1.58× end to end, 1.31× for whole
CLI cases, 1.54× for the Python API subset, and 13.29× for algorithm-only
cases.

### Carrying FASTA records across blocks

Commit `37c9828` splits a new byte block before joining the carried
bytes to its first record. It avoids copying the entire block on
standard input. Tests cover a header marker at the start of the
next block and retain the embedded-marker and CRLF cases.

Six interleaved observations per source measured 177.01 to 169.85 ms
for medium smart-gap, 281.37 to 277.35 ms for large gappy, and 335.88
to 333.27 ms for large KPIC-gappy. The small control was within 1%.
All checked output bytes matched. The change builds on the
single-byte block split from `192f57a`.

### Lowercase detection when widening packed rows

Commit `e2f8a67` checks for lowercase with the byte view used to fill the
Unicode matrix. Rows whose maximum is below ASCII `a` need no range
check. The fallback range check distinguishes lowercase from high
punctuation, and the focused tests cover that distinction alongside
NUL padding. The matrix remains writable and retains original case.

| Workload | Previous median | Candidate median |
| --- | ---: | ---: |
| Large amino-acid construction | 62.77 ms | 30.76 ms |
| Medium amino-acid construction | 31.11 ms | 18.84 ms |
| Large KPIC-gappy end to end | 332.96 ms | 288.27 ms |
| Mixed-case KPIC-smart-gap CLI | 251.04 ms | 250.41 ms |
| Dense nucleotide KPIC algorithm | 63.84 ms | 61.93 ms |

Six interleaved observations per source matched the checked outputs and
metadata exactly. Peak RSS stayed within 1%. Lowering the packed matrix
cutoff was also evaluated and slowed the small raw API control by 5.8%,
so the original packed cutoff remains.

The earlier [119-case comparison at `fef1ffa`](benchmark-results/clipkit-175x-comprehensive-fef1ffa.json)
matched `6044b1e` exactly and measured 1.67× end to end, 2.08× for the
whole CLI subset, 1.75× for Python API cases, and 13.81× for
algorithm-only cases. Tiny first-call heterotachy was its only timing
regression above 5%.

### Binary FASTA format and fallback checks

Commit `20d4a20` limits the binary block reader to streams with an
explicit FASTA format. Its text fallback starts at the original stream
position and leaves the borrowed stream open. The focused tests
exercise an explicit PHYLIP request on FASTA bytes and CRLF FASTA
starting at a nonzero stream offset. All three interleaved API cases
matched exactly and stayed within 1% in time.

A smaller-matrix bulk widening prototype left the C3, gappy, combined
API, and codon controls within 2%. The current matrix cutoff remains.

The [54-case core comparison at `fef1ffa`](benchmark-results/clipkit-175x-core-fef1ffa.json)
matched `6044b1e` exactly and measured 1.60× end to end and 13.08× for
algorithm-only cases.

### Runtime commit and rollback index

This continues the implementation history above. The [targeted samples](benchmark-results/clipkit-175x-targeted.json)
identify the reference and candidate used for each timing comparison.

| Commit | Change | Verification |
| --- | --- | --- |
| `d9385d5` | Reuse the existing matrix when all sites are kept | Exact CST output with and without logging |
| `a6edb03` | Parse common raw API inputs in memory | Exact raw and file API output, CRLF and no-tempfile tests |
| `7d08ec7` | Stream validated wide FASTA output one row at a time | Exact main and complement output and reduced peak RSS |
| `b78e9a0` | Bound integer gap checks for stop masking | Exact terminal, internal, and all-mode output and metadata |
| `f19e096` | Parse informational CLI flags before importing numerical code | Exact help/version output, import compatibility, and small CLI controls |
| `eabf605` | Widen packed ASCII to a Unicode matrix | Exact large matrix construction and trimmed output |
| `9908ef8` | Reduce histogram batch width to 256 columns | Exact sparse and dense cases with bounded RSS |
| `ae894a4` | Serialize validated FASTA rows directly | Exact headers, line wrapping, complements, and destination error handling |
| `93a6e74` | Combine excluded residue filters | Exact classification and gap summaries |
| `218368e` | Reconstruct dash counts in sparse histograms | Exact KPIC counts and dense controls |
| `c00fc29` | Use integer gap membership for medium ASCII inputs | Exact rounded fractions, ends-only results, and uppercase fallback |
| `4f44170` | Locate terminal stop positions by row | Exact masked alignments and metadata, including all-gap rows |
| `1ddf348` | Translate literal gaps during sequence-type inference | Exact type inference around the 200-residue boundary |
| `4ed4dc8` | Use byte membership for ASCII gap-only warnings | Exact log and complement, Unicode uppercase fallback |
| `c9774e8` | Share trim and log residue counts | Exact log bytes and dry-run behavior |
| `8fbeef9` | Subtract resolved counts in dash-dominated batches | Exact sparse and dense gap fractions |
| `7821042` | Build keep-only CST complements with a mask | Exact unordered and duplicate keep lists |
| `c262934` | Pass CST arrays directly into trimming | Exact million-site FASTA and unchanged public list helper |
| `9bce4f6` | Reuse the validated alignment length | Exact C3, small gappy, codon, and dense controls |
| `3a87d78` | Decode large ECOMP FASTA incrementally | Exact decoded data and trimmed output, including CRLF and Unicode headers |
| `38d8946` | Buffer per-site log file writes | Byte-identical log and alignments, handler cleanup tests |
| `a9ccb80` | Fill the final Unicode matrix row by row | Exact byte widening, padding, and truncation |
| `192f57a` | Split standard FASTA on its one-byte header marker | Exact standard and embedded-marker input |
| `913afff` | Use the block reader for files from 64 KiB | Exact small amino-acid, API, C3, and codon controls |
| `2ab1d58` | Use the block reader for raw ASCII FASTA | Exact LF and CRLF API results with no temporary input file |
| `a9e75d1` | Vectorize the existing smart-gap slope formula | Exact thresholds and duplicate-gap cases |
| `fc041aa` | Narrow histogram batches for compact small matrices | Exact small raw KPIC-smart-gap and dense controls |
| `fa671f8` | Write smaller FASTA selections directly | Exact kept and complementary bytes |
| `261b4b3` | Decode the wrapping array directly | Exact checked FASTA output |
| `fef1ffa` | Validate direct FASTA in row batches | Exact bytes and destination preservation on later invalid rows |
| `703e30d` | Reuse the FASTA wrapping buffer | Exact multiple-row output and unusual-header fallback |
| `37c9828` | Join only FASTA records spanning blocks | Exact block-boundary, embedded-marker, and CRLF input |
| `e2f8a67` | Check for lowercase on packed byte rows | Exact input case, high punctuation, NUL padding, and mixed-case controls |
| `20d4a20` | Enforce explicit binary stream formats and preserve fallback position | Exact small API controls, explicit PHYLIP rejection, and CRLF rewind |

For rollback, undo `20d4a20` before `2ab1d58` because it protects the
binary stream path introduced there. The raw FASTA cutoff also uses
`913afff`. Revert `37c9828` before `192f57a` because both change the
block-delimiter implementation. Revert `e2f8a67` before `a9ccb80`, and
`a9ccb80` before `eabf605`, to unwind the packed matrix path. Revert
`703e30d`, `fef1ffa`, `261b4b3`, and `fa671f8` in reverse history
order before removing the direct row writer in `ae894a4`. The buffered
FASTA writer in `7d08ec7` builds on the wrapping code from `684e8bb`.
Revert `8fbeef9` before `218368e` to remove the shared dash-density
sample. Revert `fc041aa` before `9908ef8` to undo the batch settings.
Revert `c262934` before `7821042` to unwind the keep-only CST path.

Set-based matching of observed states was tested on the small raw API,
small KPIC-gappy, medium nucleotide, and dense amino-acid cases. The
repeat comparison was within 2.2% for all four cases, with exact
outputs, so the existing NumPy membership check remains.

The Linux test matrix for `20d4a20` passed 552 unit and 135 integration
tests on each of CPython 3.12.14, 3.13.15, 3.14.7, and 3.15.2.
GitHub's [Version Consistency run](https://github.com/JLSteenwyk/ClipKIT/actions/runs/45685074669)
passed on branch commit `2492e7f`.

### Comprehensive checkpoint after packed case detection

The [119-case comparison at `e2f8a67`](benchmark-results/clipkit-175x-comprehensive-e2f8a67.json)
matched `6044b1e` exactly. Its unrounded end-to-end geometric mean was
approximately 1.745×. The separate first-call tiny heterotachy
case remains the only runtime regression above 5% in this checkpoint.
The raw small gappy API RSS difference was 6.9% here and zero in the
interleaved audit described above. The binary format guard and its final
suite rerun follow this checkpoint.

On Linux CPython 3.12.14, `make test` at `20d4a20` also passed
552 unit and 135 integration tests at the supported dependency floor
with NumPy 1.26.0 and Biopython 1.82.

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| `gappy` | 623.25 ms | 246.14 ms | 2.53× |
| `smart_gap` | 1,013.93 ms | 342.85 ms | 2.96× |
| `block_gappy` | 646.82 ms | 252.24 ms | 2.56× |
| `gappyout` | 1,023.80 ms | 350.09 ms | 2.92× |
| `kpi` | 755.24 ms | 326.53 ms | 2.31× |
| `kpic` | 784.87 ms | 344.26 ms | 2.28× |
| `kpi_gappy` | 646.98 ms | 291.09 ms | 2.22× |
| `kpic_gappy` | 631.18 ms | 290.44 ms | 2.17× |
| `kpi_smart_gap` | 738.11 ms | 331.28 ms | 2.23× |
| `kpic_smart_gap` | 756.53 ms | 341.52 ms | 2.22× |
| `entropy` | 1,005.06 ms | 363.07 ms | 2.77× |
| `composition_bias` | 1,035.69 ms | 363.59 ms | 2.85× |
| `heterotachy` | 1,084.93 ms | 243.67 ms | 4.45× |
| `cst` | 0.69 ms | 0.60 ms | 1.16× |
| `c3` | 6.35 ms | 5.22 ms | 1.21× |

### Major workflows

| Case | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| CLI help startup | 237.00 ms | 45.30 ms | 5.23× |
| CLI version startup | 233.53 ms | 44.79 ms | 5.21× |
| CLI small FASTA trimming | 251.82 ms | 210.58 ms | 1.20× |
| Python API from file | 4.38 ms | 2.31 ms | 1.90× |
| Python API from string | 4.71 ms | 2.29 ms | 2.06× |
| HTML report | 1,820.13 ms | 646.18 ms | 2.82× |
| Log and complementary output | 1,190.12 ms | 591.28 ms | 2.01× |
| Mixed-case amino-acid CLI | 1,050.96 ms | 247.54 ms | 4.25× |
| ECOMP input | 2.19 ms | 1.54 ms | 1.42× |
| ECOMP output | 403.95 ms | 228.57 ms | 1.77× |
| Codon trimming | 5.22 ms | 1.99 ms | 2.62× |
| Terminal stop-codon masking | 185.19 ms | 84.79 ms | 2.18× |
| Internal stop-codon masking | 181.15 ms | 84.49 ms | 2.14× |
| All stop-codon masking | 180.90 ms | 84.25 ms | 2.15× |
| One requested thread | 463.95 ms | 163.36 ms | 2.84× |
| Four requested threads | 445.84 ms | 161.53 ms | 2.76× |
| Missing ambiguity policy | 397.10 ms | 222.97 ms | 1.78× |
| Fractional ambiguity policy | 320.91 ms | 228.79 ms | 1.40× |
| Literal ambiguity policy | 320.69 ms | 228.41 ms | 1.40× |

### Geometric means

- all cases (119 cases), 4.11×
- algorithm only (47 cases), 15.25×
- all end-to-end and CLI cases (72 cases), 1.75×
- whole CLI (13 cases), 2.14×
- Python API (4 cases), 1.75×

### Import tradeoffs and C3 across input sizes

Six interleaved runs per source measured the cost of importing
`Bio.Phylo.TreeConstruction` at engine load. Eager import shortened
the timed first tiny heterotachy call from 25.64 to 4.52 ms. Whole
heterotachy CLI stayed within 1%, while the small gappy CLI rose from
211.80 to 233.31 ms. The import trace attributed about 21 ms to the
tree-construction module, and the existing heterotachy-only import
remains.

Deferring `Bio.SeqIO` after the FASTA parser and writer improvements
shortened the small gappy CLI from 217.20 to 184.37 ms. It added about
20 ms to each in-process first call for tiny FASTA and PHYLIP output.
The existing import remains. Both startup prototypes matched all
checked outputs.

| C3 case at `20d4a20` | Baseline median | Candidate median | Speedup |
| --- | ---: | ---: | ---: |
| Medium nucleotide | 6.47 ms | 5.24 ms | 1.23× |
| Codon fixture | 1.80 ms | 1.57 ms | 1.15× |
| Generated dense nucleotide (512 × 9,000) | 71.04 ms | 35.70 ms | 1.99× |

All three fixtures matched `6044b1e` exactly in six interleaved
observations per source. The larger case comes in below the 1.75×
time budget. The two smaller fixtures remain limited by input and
output work.

The [54-case core comparison at `e2f8a67`](benchmark-results/clipkit-175x-core-e2f8a67.json)
also matched the baseline exactly. It measured 1.69× end to end,
1.31× for whole CLI runs, 1.58× for Python API cases, and 14.39×
for algorithm-only cases. No core case measured more than 5% slower
in that checkpoint. The final core and comprehensive rerun at
`20d4a20` is in progress.

### Startup-bound CLI baseline audit

Six interleaved fresh-process observations per source compared the following
commands at `20d4a20` and `6044b1e`. Checked alignments matched exactly.
The last column is the reduction still needed to reach 1.75×.

| CLI workflow | Baseline median | Candidate median | Speedup | Additional reduction |
| --- | ---: | ---: | ---: | ---: |
| Small gappy | 246.78 ms | 220.06 ms | 1.12× | 79.05 ms |
| Small KPIC-gappy | 253.45 ms | 220.65 ms | 1.15× | 75.82 ms |
| Medium nucleotide smart-gap | 333.38 ms | 227.80 ms | 1.46× | 37.30 ms |
| Medium fractional ambiguity | 329.72 ms | 228.75 ms | 1.44× | 40.34 ms |
| Medium literal ambiguity | 329.02 ms | 227.34 ms | 1.45× | 39.32 ms |
| Tiny heterotachy | 243.19 ms | 238.83 ms | 1.02× | 99.86 ms |

The warmed small workloads spend only a few milliseconds trimming.
The import trace and the interleaved startup prototypes above put the
largest remaining opportunity at the Biopython import boundary. A
startup change needs to avoid both the eager tree import penalty for
ordinary CLI commands and the deferred sequence-I/O penalty for
in-process format handling. Both penalties were reproduced.

GitHub’s [Version Consistency run](https://github.com/JLSteenwyk/ClipKIT/actions/runs/45686316875)
passed on `6737aff`. The local four-version test matrix covers the
runtime code on this branch. The hosted CI workflow targets the
default branch and pull requests.

### Gather moderate-width FASTA selections as ASCII bytes

Commit `f8d239e` gathers validated residues from the native low byte
of each code point for moderate-width contiguous arrays. Signed
indices are covered in the keep and complement regression test.
Interleaved medians went from 5.24 to 4.97 ms for medium C3, 3.49 to
3.31 ms for small gappy, and 162.25 to 157.27 ms for medium smart-gap.
Checked outputs matched. The small gappy control measured 5.1% higher
peak RSS once and 3.2% lower in the repeat. C3 median RSS was unchanged.

An interleaved audit against `6044b1e` measured 1.12× for tiny FASTA
output, 2.17× for raw gappy API, and 1.26× for medium C3. All outputs
matched exactly. On Linux, `make test` at `f8d239e` passed 552 unit and
135 integration tests on each of CPython 3.12.14, 3.13.15, 3.14.7,
and 3.15.2.

### Comprehensive rerun at `20d4a20`

The [119-case comparison](benchmark-results/clipkit-175x-comprehensive-20d4a20.json)
matched `6044b1e` exactly and measured 1.757× end to end, 2.18× for
whole CLI runs, 1.78× for Python API cases, and 15.55× for algorithm-only
cases. The tiny FASTA output slowdown in that non-interleaved run
did not recur in the interleaved audit above. Tiny first-call
heterotachy remains subject to the documented import tradeoff. The
full rerun for `f8d239e` follows this checkpoint.

### Exact trimming metadata across all modes

The [metadata audit](benchmark-results/clipkit-175x-mode-metadata-f8d239e.json)
compared `f8d239e` with `6044b1e` across 23 cases spanning all 15 public
trimming modes. It checked kept and trimmed positions, rounded gap
fractions, classifications, effective gap thresholds, trimming and
stop-masking statistics, and trimmed records with IDs, names,
descriptions, and ordering. The cases include all three ambiguity
policies, codon trimming, ends-only trimming, and all three
stop-masking modes. Every check matched.

Run the script in the same Python environment as the benchmark, using
the detached baseline worktree from the reproduction commands above.

```bash
PYTHONPATH=/tmp/clipkit-175x-reference python scripts/audit_mode_metadata.py \
  /tmp/clipkit-175x-reference > /tmp/clipkit-metadata-reference.json
PYTHONPATH="$PWD" python scripts/audit_mode_metadata.py "$PWD" \
  > /tmp/clipkit-metadata-candidate.json
cmp /tmp/clipkit-metadata-reference.json /tmp/clipkit-metadata-candidate.json
```

The [54-case core rerun at `20d4a20`](benchmark-results/clipkit-175x-core-20d4a20.json)
matched exactly. Its geometric means were 1.69× end to end, 1.32× for
whole CLI runs, 1.53× for Python API cases, and 14.39× for algorithm-only
cases. No core case measured more than 5% slower. The full rerun for
the FASTA byte-gather change at `f8d239e` is in progress.

### Pushed-revision compatibility verification

On Linux x86-64 (kernel 6.8.0-55-generic), `make test` at pushed
revision `94f76c3` passed **552 unit and 135 integration tests** in
each of these environments.

| CPython | NumPy | Biopython |
| --- | --- | --- |
| 3.12.14 | 2.0.2 | 1.83 |
| 3.13.15 | 2.1.3 | 1.88 |
| 3.14.7 | 2.5.5 | 1.88 |
| 3.15.2 | 2.5.5 | 1.88 |
| 3.12.14 | 1.26.0 | 1.82 |

The final row uses the supported dependency floors. The 23-case exact
metadata audit also matched the baseline in that environment. The
[GitHub Version Consistency run](https://github.com/JLSteenwyk/ClipKIT/actions/runs/45687832853)
passed on `94f76c3`. The macOS test workflow requires a default-branch
push or pull request. No hosted macOS results are available for this
branch.

The benchmark worker uses `f8d239e`, whose production code, tests, and
benchmark runner are identical to those at `94f76c3`.

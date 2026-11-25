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

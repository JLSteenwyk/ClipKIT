# Changelog

## Unreleased

- Added `--ambiguity_handling {missing,fractional,literal}` (and the equivalent Python API argument) for explicit IUPAC ambiguity handling.
- Made conservative `missing` handling the default: ambiguity symbols are excluded from entropy, composition, heterotachy, and KPI/KPIC state counts, and contribute to the effective unavailable fraction in gap-based modes.
- Added optional fractional IUPAC weighting for entropy, composition, and heterotachy's clade-entropy calculation while keeping KPI/KPIC classification conservative, plus a `literal` compatibility mode for the legacy ambiguity interpretation.
- Improved sequence-type auto-detection for IUPAC-rich nucleotide alignments and added separate configured-gap, ambiguity, and resolved-state tracks to HTML reports.
- Documented nucleotide/protein ambiguity maps, configured-gap precedence, mode-specific behavior, and unchanged alignment output.

## 2.13.2

- Expanded behavioral test coverage from 84.53% to 98.83% line coverage, increasing the suite from 552 to 619 tests without modifying production behavior.
- Added coverage for public Python APIs, CLI validation and execution paths, file output, report generation, malformed inputs, and reusable in-process execution.
- Added boundary and regression tests for alignment parsing and trimming, codon and stop-codon handling, ECOMP archives, smart-gap selection, and gappyout behavior.
- Exercised benchmark workers, fresh-process execution, CLI/API equivalence, suite orchestration, result comparison, and report generation.

## 2.13.1

- Accelerated terminal, internal, and all stop-codon masking by avoiding redundant full-matrix Unicode normalization for inputs already known to be uppercase.
- Vectorized entropy and composition-bias scoring across alignment columns while preserving exact four-decimal scores and trimming decisions.
- Vectorized codon-site expansion for faster codon-aware trimming on alignments with many selected sites.
- Expanded reproducible benchmarks to cover every trimming mode, CLI/API paths, threads, stop-codon modes, supported formats, exact baseline comparison, warm-ups, runtime ranges, and peak memory.
- Accelerated gappy, smart-gap, KPI, and KPIC column counting by combining and fusing temporary-array encoding passes without changing site statistics or trimming decisions.
- Reused exact classification counts for uppercase KPI/KPIC gap-combination modes, avoiding a second full alignment scan while retaining separate raw and normalized passes for mixed-case inputs.
- Added a focused 54-case core-mode benchmark matrix with dense/sparse AA and NT fixtures, exact statistic/position hashes, resolved thresholds, effective thread counts, and retained JSON/CSV reports.

## 2.13.0

- Added `--remove_stop_codons {terminal,internal,all}` for masking selected in-frame DNA or RNA stop codons as gaps before alignment statistics and codon-aware trimming.
- Added the same stop-codon masking modes to the public Python API, with validation for nucleotide input, codon mode, and triplet alignment lengths.
- Added human-readable and JSON masking summaries with separate terminal and internal stop-codon counts.

## 2.12.2

- Replaced per-column character sorting with batched compact counting, substantially accelerating KPI/KPIC, entropy, composition-bias, and gap-based trimming while preserving exact trimming decisions.
- Reused cached per-site gap statistics during smart-gap threshold selection instead of scanning the alignment twice.
- Accelerated alignment output by materializing complete NumPy rows directly as sequence strings instead of allocating one Python object per character.
- Expanded the benchmark suite to cover representative small, medium, and large workloads with output hashes and peak-memory measurements.

## 2.12.1

- Improved MSA construction from BioPython alignments by converting sequence strings directly into a NumPy character matrix instead of first building per-character Python lists.
- Bumped the package and Galaxy wrapper version to `2.12.1` because `2.12.0` already exists on PyPI and PyPI release files are immutable.
- Updated manual release instructions to use `python -m build` and `twine check` before upload.
- Removed the stale `bdist_wheel --universal` release command so ClipKIT publishes Python 3-only wheels for its supported Python 3.10+ range.
- Excluded test modules from installable wheel packages.

## 2.12.0

- When `-g`/`--gaps` is provided without `-m`/`--mode`, the default mode is now `gappy` instead of `smart-gap`. This ensures user-specified gap thresholds are honoured rather than silently overridden by dynamic threshold calculation.
- A warning is now emitted when `-g`/`--gaps` is explicitly combined with a mode that ignores it (`smart-gap`, `kpi-smart-gap`, `kpic-smart-gap`, or `gappyout`).
- Updated documentation to reflect the new `-g`/`-m` interaction behaviour.

## 2.11.4

- Added performance-focused thread auto-tuning for KPI/KPIC-family modes to reduce overhead on workloads where fewer threads are faster.
- Updated CLI help text to clarify that `--threads` is a requested thread count and may be tuned downward in KPI/KPIC-family modes.
- Added `--dry_run` and `--validate_only` CLI modes for safer preview/validation workflows.
- Added `--report_json` to emit machine-readable run summaries (with optional default path behavior).
- Added `--plot_trim_report` to generate an interactive HTML report with per-site tracks and trimmed-column highlighting.
- Added export controls in the trim report to save per-site tracks and alignment preview as PNG images.
- Added a new `gappyout` trimming mode with automatic, gap-distribution-based threshold selection (gappyout-inspired behavior).
- Added a new `block-gappy` trimming mode for trimming contiguous runs of high-gappyness sites.
- Added a new `composition-bias` trimming mode for trimming sites with strong compositional skew.
- Added a new `heterotachy` trimming mode that infers a parsimony guide tree and trims sites with high clade-level entropy variation.
- Dropped Python 3.9 support and set the minimum supported version to Python 3.10.
- Improved runtime performance by reducing duplicate alignment matrix construction and optimizing frequency/classification hot paths.
- Hardened CLI/API execution behavior and logging lifecycle for repeated in-process runs.
- Added regression tests for complementary output handling, thread validation, invalid input handling, logger/handler cleanup, and thread heuristic behavior.
- Switched release workflow to token-based PyPI publishing (`PYPI_API_TOKEN`) to match the project's standard release flow and avoid trusted-publisher configuration failures.

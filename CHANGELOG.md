# Changelog

## 2.11.1

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
- Fixed release workflow trusted-publishing environment claim (`environment: pypi`) for PyPI publishing.

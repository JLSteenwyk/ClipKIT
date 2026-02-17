# Changelog

## 2.10.3

- Added performance-focused thread auto-tuning for KPI/KPIC-family modes to reduce overhead on workloads where fewer threads are faster.
- Updated CLI help text to clarify that `--threads` is a requested thread count and may be tuned downward in KPI/KPIC-family modes.
- Added `--dry_run` and `--validate_only` CLI modes for safer preview/validation workflows.
- Added `--report_json` to emit machine-readable run summaries (with optional default path behavior).
- Added a new `gappyout` trimming mode with automatic, gap-distribution-based threshold selection (gappyout-inspired behavior).
- Dropped Python 3.9 support and set the minimum supported version to Python 3.10.
- Improved runtime performance by reducing duplicate alignment matrix construction and optimizing frequency/classification hot paths.
- Hardened CLI/API execution behavior and logging lifecycle for repeated in-process runs.
- Added regression tests for complementary output handling, thread validation, invalid input handling, logger/handler cleanup, and thread heuristic behavior.

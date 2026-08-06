.. _change_log:


Change log
==========

Major changes to ClipKIT are summarized here.

**Unreleased**

**2.14.0**

- Added ``--ambiguity_handling {missing,fractional,literal}`` and equivalent
  Python API support for explicit IUPAC ambiguity handling.
- Made conservative ``missing`` handling the default: ambiguity symbols are
  excluded from entropy, composition, heterotachy, and KPI/KPIC state counts,
  and contribute to the effective unavailable fraction in gap-based modes.
- Added fractional IUPAC weighting for entropy, composition, and heterotachy's
  clade-entropy calculation while
  keeping KPI/KPIC classification conservative, plus a ``literal`` legacy
  compatibility mode.
- Improved auto-detection of IUPAC-rich nucleotide alignments and added
  configured-gap, ambiguity, and resolved-state tracks to HTML reports.
- Documented the supported nucleotide/protein ambiguity maps, configured-gap
  precedence, mode-specific behavior, and unchanged alignment output.

**2.13.2**

- Expanded behavioral test coverage from 84.53% to 98.83% line coverage,
  increasing the suite from 552 to 619 tests without modifying production
  behavior.
- Added coverage for public Python APIs, CLI validation and execution paths,
  file output, report generation, malformed inputs, and reusable in-process
  execution.
- Added boundary and regression tests for alignment parsing and trimming,
  codon and stop-codon handling, ECOMP archives, smart-gap selection, and
  gappyout behavior.
- Exercised benchmark workers, fresh-process execution, CLI/API equivalence,
  suite orchestration, result comparison, and report generation.

**2.13.1**

- Accelerated terminal, internal, and all stop-codon masking by skipping
  redundant full-matrix Unicode normalization for inputs already known to be
  uppercase.
- Vectorized entropy and composition-bias scoring across alignment columns
  while preserving exact four-decimal scores and trimming decisions.
- Vectorized codon-site expansion for faster codon-aware trimming when many
  sites are selected.
- Expanded reproducible benchmarks across every trimming mode, CLI/API paths,
  thread counts, stop-codon modes, and supported formats, with exact baseline
  comparison, warm-ups, runtime ranges, and peak-memory reporting.
- Accelerated gappy, smart-gap, KPI, and KPIC column counting by combining and
  fusing temporary-array encoding passes without changing site statistics or
  trimming decisions.
- Reused exact classification counts for uppercase KPI/KPIC gap-combination
  modes while retaining separate raw and normalized passes for mixed-case
  inputs.
- Added a focused 54-case core-mode benchmark matrix spanning dense and sparse
  amino-acid and nucleotide inputs, exact result hashes, resolved thresholds,
  effective thread counts, and retained JSON/CSV reports.

**2.13.0**

- Added ``--remove_stop_codons {terminal,internal,all}`` to mask selected
  in-frame DNA or RNA stop codons before gap statistics and codon-aware
  trimming.
- Added equivalent public Python API support and validation for nucleotide,
  codon-mode, and triplet-length requirements.
- Added terminal/internal masking counts to human-readable and JSON reports.

**2.12.2**
Improved trimming and output performance while preserving exact results:

- Replaced per-column character sorting with batched compact counting for
  faster KPI/KPIC, entropy, composition-bias, and gap-based trimming.
- Reused cached per-site gap statistics during smart-gap threshold selection.
- Materialized complete NumPy rows directly as sequence strings during output
  instead of allocating one Python object per character.
- Expanded benchmarks to cover small, medium, and large workloads with output
  hashes and peak-memory measurements.

**2.12.1**
Improved performance and packaging release behavior:

- Faster MSA construction from BioPython alignments by converting sequence
  strings directly into a NumPy character matrix instead of first building
  per-character Python lists.
- Bumped the package and Galaxy wrapper version to ``2.12.1`` because
  ``2.12.0`` already exists on PyPI and PyPI release files are immutable.
- Updated manual release instructions to use ``python -m build`` and
  ``twine check`` before upload.
- Removed the stale ``bdist_wheel --universal`` release command so ClipKIT
  publishes Python 3-only wheels for its supported Python 3.10+ range.
- Excluded test modules from installable wheel packages.

**2.12.0**
Improved ``-g``/``-m`` argument interaction:

- When ``-g``/``--gaps`` is provided without ``-m``/``--mode``, the default mode is now ``gappy`` instead of ``smart-gap``, so user-specified gap thresholds are honoured rather than silently overridden.
- A warning is now emitted when ``-g``/``--gaps`` is explicitly combined with a mode that ignores it (``smart-gap``, ``kpi-smart-gap``, ``kpic-smart-gap``, or ``gappyout``).
- Updated documentation to reflect the new ``-g``/``-m`` interaction behaviour.

**2.11.4**
Improved performance and usability through:
- KPI/KPIC-family thread auto-tuning to avoid slowdowns from excessive parallelism
- Additional performance optimizations in hot alignment-processing paths
- Updated CLI/documentation text clarifying `--threads` behavior as a requested value
- Added CLI safety/reporting options: ``--dry_run``, ``--validate_only``, and ``--report_json``
- Added ``--plot_trim_report`` for interactive HTML visualization of per-site metrics and trimmed columns
- Added trim-report export controls for saving per-site tracks and alignment preview as PNG files
- Added a ``gappyout`` trimming mode with automatic, gap-distribution-based threshold selection
- Added a ``block-gappy`` trimming mode for contiguous runs of high-gappyness sites
- Added a ``composition-bias`` trimming mode for sites with strong compositional skew
- Added a ``heterotachy`` trimming mode using a parsimony guide tree and clade-level entropy variation
- Dropped Python 3.9 support and raised the minimum supported version to Python 3.10
- Robustness improvements for repeated in-process execution and logging behavior
- Switched release workflow to token-based PyPI publishing (``PYPI_API_TOKEN``) to match the project's standard release flow and avoid trusted-publisher configuration failures

**2.7.0**
Significant performance improvements with 2-60x speedup (average 9x) through:
- Vectorized NumPy operations for column frequency calculations
- Parallel processing with adaptive thread allocation
- Optimized memory usage with batch processing
- Caching for expensive computations
- Added support for Python 3.12 and 3.13

**2.4.1**
The ends_only function now handles the case where nothing gets trimmed.

**2.4.0**
Added a new function called ends_only, which removes sites that would
be trimmed for a given mode, but only at the ends of the alignment.
For example, if the sites that should be trimmed include
[0, 1, 2, 4, 5, 6, 14, 15, 16] for smart-gap mode and an alignment of
length 16, adding the ends_only mode will result in [0, 1, 2, 14, 15, 16]
being the sites that will be trimmed. Specify this argument with -eo, \-\-ends_only.

**2.3.0**
Added support for Python version 3.11

**2.2.3**
Fixed gap character handling. The help message was incongruent
with what was happening underneath the hood.

**2.1.2**
Incorporate codon-based trimming. When one position in a codon gets trimmed based on the mode
being used, the whole codon will get trimmed from the alignment.

**1.4.0**
new argument for specifying if sequences are amino acids or nucleotides

**1.3.0**
long description of sequences, rather than identifiers, are kept in the ClipKIT output

**1.1.5**
carried over code base to biopython, v1.79

**1.1.0:**
smart-gap trimming is introduced and is now the default trimming approach used in ClipKIT.
smart-gap trimming is a dynamic approach to determine the appropriate gaps threshold for an alignment.

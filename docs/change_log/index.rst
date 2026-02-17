.. _change_log:


Change log
==========

Major changes to ClipKIT are summarized here.

**2.10.0**
Improved performance and usability through:
- KPI/KPIC-family thread auto-tuning to avoid slowdowns from excessive parallelism
- Additional performance optimizations in hot alignment-processing paths
- Updated CLI/documentation text clarifying `--threads` behavior as a requested value
- Added CLI safety/reporting options: ``--dry_run``, ``--validate_only``, and ``--report_json``
- Robustness improvements for repeated in-process execution and logging behavior

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

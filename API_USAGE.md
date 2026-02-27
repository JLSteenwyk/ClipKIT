# ClipKIT Python API Usage

ClipKIT can be used programmatically within Python scripts. This document describes how to use the ClipKIT API.

## Installation

First, install ClipKIT:
```bash
pip install clipkit
```

## Basic Usage

```python
from clipkit import clipkit
from clipkit.modes import TrimmingMode
from clipkit.files import FileFormat

# Trim an alignment file
trim_run, stats = clipkit(
    input_file_path="alignment.fasta",
    output_file_path="trimmed_alignment.fasta",
    mode=TrimmingMode.smart_gap,
)

# Access trimming statistics
print(f"Sites kept: {stats.output_length}")
print(f"Sites trimmed: {stats.trimmed_length}")
```

## Parameters

### Required Parameters
- `input_file_path` or `raw_alignment`: Input alignment (file path or string)

### Optional Parameters
- `output_file_path`: Output file path (if not specified, returns `trim_run` + `stats`)
- `mode`: Trimming mode (default: `TrimmingMode.smart_gap`)
  - Options: `smart_gap`, `entropy`, `gappy`, `block_gappy`, `gappyout`, `composition_bias`, `kpic`, `kpic_smart_gap`, `kpic_gappy`, `kpi`, `kpi_smart_gap`, `kpi_gappy`, `cst`, `c3`
- `gaps`: Threshold in `[0,1]` (default: `0.9`; auto-calculated for `smart_gap` and `gappyout`; default `0.8` for `entropy` and `composition_bias`)
  - Interpreted as gappyness for gap-based modes, normalized Shannon entropy for `entropy`, and normalized compositional skew for `composition_bias`
- `gap_characters`: List of gap characters (default: auto-detect based on sequence type)
- `input_file_format`: Input format (default: `FileFormat.fasta`)
- `output_file_format`: Output format (default: `FileFormat.fasta` for API usage)
- `sequence_type`: Sequence type - `SeqType.aa` or `SeqType.nt` (default: auto-detect)
- `codon`: Enable codon-based trimming (default: `False`)
- `ends_only`: Trim only alignment ends (default: `False`)
- `threads`: Requested number of threads for parallel processing (default: `1`)
  - **Note**: ClipKIT may use fewer threads for KPI/KPIC-family modes when that is expected to be faster
- `plot_trim_report_path`: Optional HTML report path with per-site tracks and trimmed-site highlighting

## Examples

### Parallel Processing for Large Alignments
```python
from clipkit import clipkit

# Use multiple threads for large alignments
trim_run, stats = clipkit(
    input_file_path="large_alignment.fasta",
    output_file_path="trimmed_large.fasta",
    threads=8  # Use 8 threads for faster processing
)
```

### Different Trimming Modes
```python
from clipkit import clipkit
from clipkit.modes import TrimmingMode

# Keep only parsimony-informative sites
trim_run, stats = clipkit(
    input_file_path="alignment.fasta",
    mode=TrimmingMode.kpi,
)

# Gappy mode with custom threshold
trim_run, stats = clipkit(
    input_file_path="alignment.fasta",
    mode=TrimmingMode.gappy,
    gaps=0.5,  # Remove sites with >50% gaps
)

# Gappyout-inspired mode with automatic threshold selection
trim_run, stats = clipkit(
    input_file_path="alignment.fasta",
    mode=TrimmingMode.gappyout,
)

# Write an interactive trim report
trim_run, stats = clipkit(
    input_file_path="alignment.fasta",
    mode=TrimmingMode.gappy,
    gaps=0.5,
    plot_trim_report_path="alignment.trim_report.html",
)
```

### Working with Raw Alignment Strings
```python
from clipkit import clipkit

alignment_string = ">seq1\nATGC--\n>seq2\nATG---\n>seq3\nATGCAT"

trim_run, stats = clipkit(
    raw_alignment=alignment_string,
    output_file_path="output.fa",
)
```

## Performance Considerations

- For alignments with <5000 sites, single-threaded mode is typically faster
- For large alignments, using multiple threads can provide significant speedup
- The optimal number of threads depends on your system and alignment size
- For KPI/KPIC-family modes, ClipKIT may automatically use fewer threads than requested to reduce overhead
- Results are identical regardless of thread count (fully reproducible)

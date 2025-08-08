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
print(f"Sites kept: {stats.sites_kept}")
print(f"Sites trimmed: {stats.sites_trimmed}")
```

## Parameters

### Required Parameters
- `input_file_path` or `raw_alignment`: Input alignment (file path or string)

### Optional Parameters
- `output_file_path`: Output file path (if not specified, returns Bio.Align object)
- `mode`: Trimming mode (default: `TrimmingMode.smart_gap`)
  - Options: `smart_gap`, `gappy`, `kpic`, `kpic_smart_gap`, `kpic_gappy`, `kpi`, `kpi_smart_gap`, `kpi_gappy`, `c3`
- `gaps`: Gap threshold (default: 0.9 for gappy mode, auto-calculated for smart-gap)
- `gap_characters`: List of gap characters (default: auto-detect based on sequence type)
- `input_file_format`: Input format (default: `FileFormat.fasta`)
- `output_file_format`: Output format (default: same as input)
- `sequence_type`: Sequence type - `SeqType.aa` or `SeqType.nt` (default: auto-detect)
- `codon`: Enable codon-based trimming (default: `False`)
- `ends_only`: Trim only alignment ends (default: `False`)
- `threads`: Number of threads for parallel processing (default: `1`)
  - **Note**: Parallel processing is only activated for alignments >5000 sites

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
- For large alignments (>5000 sites), using multiple threads can provide significant speedup
- The optimal number of threads depends on your system and alignment size
- Results are identical regardless of thread count (fully reproducible)
from dataclasses import dataclass

@dataclass
class TrimmingStats:
    alignment_length: int
    output_length: int
    trimmed_length: int
    trimmed_percentage: float

def calculate_stats(alignment, keepD, trimD):
    alignment_length = alignment.get_alignment_length()
    output_length = len(next(iter(keepD.values())))
    trimmed_length = len(next(iter(trimD.values())))
    stats = {
        "alignment_length": alignment_length,
        "output_length": output_length,
        "trimmed_length": trimmed_length,
        "trimmed_percentage": round((trimmed_length / alignment_length) * 100, 3),
    }
    return TrimmingStats(**stats)

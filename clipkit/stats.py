def calculate_stats(alignment, keepD, trimD):
    alignment_length = alignment.get_alignment_length()
    output_length = len(next(iter(keepD.values())))
    trimmed_length = len(next(iter(trimD.values())))
    return {
        "alignment_length": alignment_length,
        "output_length": output_length,
        "trimmed_length": trimmed_length,
        "trimmed_percentage": round((trimmed_length / alignment_length) * 100, 3),
    }

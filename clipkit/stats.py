from dataclasses import dataclass

from Bio.Align import MultipleSeqAlignment

@dataclass
class TrimmingStats:
    alignment: MultipleSeqAlignment
    keepD: dict
    trimD: dict

    def alignment_length(self):
        return self.alignment.get_alignment_length()

    def output_length(self):
        return len(next(iter(self.keepD.values())))

    def trimmed_length(self):
        return len(next(iter(self.trimD.values())))

    def trimmed_percentage(self):
        return round((self.trimmed_length / self.alignment_length) * 100, 3)

    def summary(self):
        return {
            "alignment_length": self.alignment_length,
            "output_length": self.output_length,
            "trimmed_length": self.trimmed_length,
            "trimmed_percentage": self.trimmed_percentage,
        }

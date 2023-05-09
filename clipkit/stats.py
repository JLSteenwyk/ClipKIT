from dataclasses import dataclass

from Bio.Align import MultipleSeqAlignment

from .msa import MSA


@dataclass
class TrimmingStats:
    alignment: MultipleSeqAlignment
    keepMSA: MSA
    trimMSA: MSA

    @property
    def alignment_length(self):
        return self.alignment.get_alignment_length()

    @property
    def output_length(self):
        return self.keepMSA.length

    @property
    def trimmed_length(self):
        return self.alignment_length - self.output_length

    @property
    def trimmed_percentage(self):
        return round((self.trimmed_length / self.alignment_length) * 100, 3)

    @property
    def summary(self):
        return {
            "alignment_length": self.alignment_length,
            "output_length": self.output_length,
            "trimmed_length": self.trimmed_length,
            "trimmed_percentage": self.trimmed_percentage,
        }

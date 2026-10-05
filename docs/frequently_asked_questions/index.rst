.. _faq:


FAQ
===

**Does ClipKIT trim amino acids, nucleotides, or codons?**

ClipKIT trims amino acid and nucleotide alignments. For codon-aligned nucleotide
MSAs, ``--codon`` keeps trimming decisions in whole codons. The optional
``--remove_stop_codons`` mode can first mask terminal, internal, or all in-frame
stop codons as gaps.

|

**Which trimming mode should I use before selection analyses (e.g., HyPhy or PAML)?**

For codon-aligned nucleotide MSAs used to estimate dN/dS or test for selection,
we recommend a gap-based mode with codon-aware trimming, such as
``clipkit coding.fa --codon --sequence_type nt --remove_stop_codons all`` (the
default ``smart-gap`` mode) or ``-m gappy -g <threshold>`` for explicit control
over gap tolerance. Lower ``-g`` values trim more aggressively. We advise against
``kpi`` and ``kpic`` modes for these analyses: singleton sites are real
substitutions, often on terminal branches, and are informative for selection
inference. In codon mode, ``kpic`` modes also remove every codon that contains a
singleton site, including its constant and parsimony-informative positions. See
`Codon`_ for details.

.. _`Codon`: ../advanced/index.html#codon

|

**Is there a website application of ClipKIT?**

Currently, ClipKIT is only a command line tool.

|

**If tree inference with no trim works well, why even trim?**

Tree inference with trimmed multiple sequence alignments is computationally efficient.
In other words, shorter alignments require less computational time and memory during tree
search. We found that ClipKIT reduced computation time by an average of 20%. As datasets
continuously become bigger, an alignment trimming algorithm that can reduce computational
time will be of great value. 

|

**What characters are considered gaps?**

For amino acids ?, \*, -, X; for nucleotides, the same characters and N.

|

**How are ambiguous bases or amino acids handled?**

By default, recognized IUPAC ambiguity symbols are treated as missing evidence:
they are excluded from entropy/composition and KPI/KPIC state counts and are
included in the effective unavailable fraction used by gap-based modes. Use
``--ambiguity_handling fractional`` to distribute them equally among possible
states for entropy/composition calculations, or ``literal`` for the legacy
ambiguity interpretation. Configured gap characters take precedence, and
ClipKIT does not rewrite the input symbols. See `Ambiguity handling`_ for the
complete mappings and mode-specific details.

.. _`Ambiguity handling`: ../advanced/index.html#ambiguity-handling

|

**I am having trouble install ClipKIT, what should I do?**

Please install ClipKIT using a virtual environment as directed in the installation instructions.
If you are still running into issues after installing in a virtual environment, please contact the
main software developer via email_ or twitter_.

.. _email: https://jlsteenwyk.com/contact.html
.. _twitter: https://twitter.com/jlsteenwyk

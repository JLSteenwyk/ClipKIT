.. _change_log:


Change log
==========

^^^^^

Major changes to ClipKIT are summarized here.

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

.. image:: ../_static/img/ClipKIT_logo_top_only_v1.jpg
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/ClipKIT

.. _faq:


FAQ
===

**Does ClipKIT trim amino acids, nucleotides, or codons?**

ClipKIT trims amino acid and nucleotide alignments. Currently, ClipKIT does not trim codons. 

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

^^^^^
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

|

**What characters are considered gaps?**

For amino acids ?, \*, -, X; for nucleotides, the same characters and N.

|

**I am having trouble install ClipKIT, what should I do?**

Please install ClipKIT using a virtual environment as directed in the installation instructions.
If you are still running into issues after installing in a virtual environment, please contact the
main software developer via email_ or twitter_.

.. _email: https://jlsteenwyk.com/contact.html
.. _twitter: https://twitter.com/jlsteenwyk

^^^^^

Performance Assessment
======================

^^^^^


In brief, performance assessment and comparison of multiple trimming alignment software
revealed that ClipKIT with nearly any mode and trimAl with the 'gappyout' parameter are
the top-performing software. Here, we provide greater detail into the empirical datasets
used to assess alignment trimming performance.

.. image:: ../_static/img/Performance_summary_desirability.jpg

**Summary Figure. Rank-based assessment of ClipKIT's performance compared to other 
alignment trimming software revealed ClipKIT is a top-performing software.** A dataset of
2,002 amino acid (AA) and nucleotide (NT) alignments across 24 Mammals (A, B) and 2,832
AA and NT alignments across 12 yeasts (C, D) were trimmed using ClipKIT,
`trimAl <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2712344/>`_, 
`BMGE <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3017758/>`_,
`Noisy <https://www.ncbi.nlm.nih.gov/pubmed/18577231>`_, and
`Gblocks <https://www.ncbi.nlm.nih.gov/pubmed/17654362>`_. The resulting trimmed 
alignments were then used to infer the evolutionary history of each alignment using
`IQ-Tree <https://www.ncbi.nlm.nih.gov/pubmed/32011700>`_. To assess the accuracy and
support of each topology, we calculated normalized Robinson-Foulds distances (nRF) and
average bootstrap support (ABS). For every gene, we obtained one performance score by
taking the average of the desirability transformed nRF and ABS values. Using the
resulting desirability values, trimmers were ranked from 1 to 9 where lower ranks
indicate better performance. We then determined the best performing algorithm per data
set and plotted them in order of best to worst. We found that ClipKIT, trimAl (with
the 'gappyout' parameter), and no trimming performed the best across all data matrices.
Abbreviations of trimmers and parameters are as follows: ClipKIT: kg = kpi-gappy mode;
ClipKIT: k = kpi mode; ClipKIT: g = gappy mode; trimAl: go = gappyout; trimAl: s =
strictplus.

|

^^^^^

.. image:: ../_static/img/Alignment_length_zscores.jpg

**Figure 1. ClipKIT provides 'light' and 'heavy' trimming options that perform
well.** Using the same dataset of sequence alignments for mammals (A, B) and yeast (C, D),
comparisons of alignment lengths, we calculated Z-scores of alignment lengths for every 
gene. ClipKIT's 'gappy' mode conducts the lightest trimming followed by trimAl with the 
'gappyout' parameter. Interestingly, ClipKIT's 'kpi-gappy' and 'kpi' modes trim a 
significant amount of the alignment compared to other trimmers. Despite alignment length 
being associated with strong phylogenetic signal, ClipKIT's 'kpi-gappy' and 'kpi' modes 
resulted in alignments that were still rich in phylogenetic information (see Summary 
Figure.) The red dashed line represents a Z-score of zero. Abbreviations of trimmers and 
parameters is the same as the summary figure.

|

^^^^^

.. image:: ../_static/img/RF_and_ABS_zscores.jpg

**Figure 2. ClipKIT is a top-performing software across metrics associated with strong
phylogenetic signal.** Using the same dataset of sequence alignments for mammals and
yeast, nRF distances (A-D) and ABS (E-H; using IQ-Tree's UFBoot), were calculated on a
per-gene basis. For each gene, nRF and ABS values were z-transformed. Across these metrics
ClipKIT -- no matter the mode -- was a top-performing software.
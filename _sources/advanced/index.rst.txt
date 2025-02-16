Advanced Usage
==============

^^^^^

This section describes the various features and options of ClipKIT.

- Modes_
- Output_
- Log_
- Complementary_
- Codon_
- `Custom site trimming (cst mode)`_
- Gaps_
- `Gap Characters`_
- `Sequence Type`_
- `All options`_

|

.. _Modes:

Modes
-----

Herein, we describe the various trimming modes implemented in ClipKIT. If you are unsure which is appropriate for you,
**we recommend using the default smart-gap trimming mode**. 

ClipKIT can be run with eight different modes, which are specified with the -m/--mode argument.
*Default: 'smart-gap'*

* smart-gap: dynamic determination of gaps threshold
* gappy: trim all sites that are above a threshold of gappyness (default: 0.9)
* kpic: keep only parismony informative and constant sites
* kpic-smart-gap: a combination of kpic- and smart-gap-based trimming 
* kpic-gappy: a combination of kpic- and gappy-based trimming
* kpi: keep only parsimony informative sites
* kpi-smart-gap: a combination of kpi- and smart-gap-based trimming
* kpi-gappy: a combination of kpi- and gappy-based trimming
* c3: remove third codon position from alignment
* cst: custom site trimming (remove sites specified by the user)

.. code-block:: shell

	# smart-gap-based trimming
	clipkit <input>
	clipkit -m smart-gap

	# gappy-based trimming
	clipkit <input> -m gappy

	# kpic-based trimming
	clipkit <input> -m kpic

	# kpic- and smart-gap-based trimming
	clipkit <input> -m kpic-smart-gap

	# kpic- and gappy-based trimming
	clipkit <input> -m kpic-gappy

	# kpi-based trimming
	clipkit <input> -m kpi

	# kpi- and smart-gap-based trimming
	clipkit <input> -m kpi-smart-gap

	# kpi- and gappy-based trimming
	clipkit <input> -m kpi-gappy

	# remove third codon position
	clipkit <input> -m c3

	# conduct site-specific trimming
	clipkit <input> -m cst -a <auxiliary file>

.. _Output:

|

Output
------

By default, output files will have the same name as the input file with the suffix ".clipkit"
appended to the name. Users can specify output file names with the -o option. 

.. code-block:: shell

	# specify output
	clipkit <input> -o <output>

|

.. _Log:

Log
---
It can be very useful to have information about the each position in an alignment. For
example, this information could be used in alignment diagnostics, fine-tuning of trimming
parameters, etc. To create the log file, use the -l/\\-\\-log option. Using this option
will create a four column file with the suffix 'clipkit.log'. *Default: off*

* col1: position in the alignment (starting at 1)
* col2: reports if site was trimmed or kept (trim or keep, respectively)
* col3: reports if the site is parsimony informative or not (PI or nPI, respectively)
* col4: reports the gappyness of the position (number of gaps / entries in alignment)

.. code-block:: shell

	clipkit <input> -l 

|

.. _Complementary:

Complementary
-------------

Having an alignment of the sequences that were trimmed can be useful for other analyses. 
To obtain an alignment of the sequences that were trimmed, use the -c/\\-\\-complementary 
option.

.. code-block:: shell

	clipkit <input> -c

Output file with the suffix '.clipkit.complementary'

|

.. _Codon:

Codon
-----

Trims codon-based alignments. If one position in a codon should be trimmed, the whole
codon will be trimmed. To conduct codon-based trimming, use the -co/\\-\\-codon argument.

.. code-block:: shell

	clipkit <input> --codon

    # or

	clipkit <input> --co

|


.. _`Custom site trimming (cst mode)`:

Custom site trimming (cst mode)
-------------------------------

Custom site trimming specified using a tab-delimited text file specified using the -a argument.

.. code-block:: shell

	clipkit <input> -m cst -a <auxiliary_file>

|

The `auxiliary_file` is a two column tab-delimited file wherein the first column is the site
(starting at 1) and the second column specifies if the site should be kept or trimmed using the
strings "keep" or "trim".

.. code-block:: shell

	cat auxiliary_file.txt

	1	keep
	2	trim
	3	keep
	4	keep
	5	keep
	6	keep

|

Alternatively, users can specify sites that are only kept or trimmed using the `auxiliary_file`.
For example, the following would be equivalent to the auxiliary file described above.

.. code-block:: shell

	cat auxiliary_file.txt

	2	trim

|

Similarly, the following would conduct the trimming, wherein the second site is removed but all
others are kept. 

.. code-block:: shell

	cat auxiliary_file.txt

	1	keep
	3	keep
	4	keep
	5	keep
	6	keep

|

.. _Gaps:

Gaps
----

Positions with gappyness greater than threshold will be trimmed. 
Must be between 0 and 1. (Default: 0.9). This argument is ignored
when using the kpi and kpic mdoes of trimming as well as an 
iteration of trimming that uses smart-gap.

To specify a gaps threshold, use the -g/\\-\\-gaps argument.

.. code-block:: shell

	clipkit <input> --gaps 0.4

    # or

	clipkit <input> --g 0.4

|

.. _`Gap Characters`:

Gap Characters
--------------

Specifies gap characters used in the input file. For example,
"NnXx-?" would specify that "N", "n", "X", "x", "-", and "?" are
gap characters. Note, the first gap character cannot be "-" because
the parser will interpret the gaps list as a new argument.

.. code-block:: shell

	clipkit <input> -gc NnXx-?

|

.. _`Sequence Type`:

Sequence Type
-------------

Specifies the type of sequences in the input file. The default
is auto-detection of sequence type. Valid options
include aa or nt for amino acids and nucleotides. This argument
is case insensitive. This matters for what characters are
considered gaps. For amino acids, -, ?, \*, and X are considered
gaps. For nucleotide sequences, the same characters are
considered gaps as well as N.

.. code-block:: shell

	clipkit <input> -s aa

Specify input sequences are amino acids

.. code-block:: shell

	clipkit <input> -s nt

Specify input sequences are nucleotides 

|

.. _`All options`:

All options
---------------------


+-----------------------------+-------------------------------------------------------------------------+
| Option                      | Usage and meaning                                                       |
+=============================+=========================================================================+
| -h/\\-\\-help               | Print help message                                                      |
+-----------------------------+-------------------------------------------------------------------------+
| -v/\\-\\-version            | Print software version                                                  |
+-----------------------------+-------------------------------------------------------------------------+
| -m/\\-\\-mode               | Specify trimming mode (default: smart-gap)                              |
+-----------------------------+-------------------------------------------------------------------------+
| -o/\\-\\-output             | Specify output file name                                                |
+-----------------------------+-------------------------------------------------------------------------+
| -g/\\-\\-gaps               | Specify gappyness threshold (between 0 and 1). *Default: 0.9*           |
+-----------------------------+-------------------------------------------------------------------------+
| -gc/\\-\\-gap_characters    | Specifies gap characters used in input file (AAs: Xx-?*; NTs: XxNn-?*   |
+-----------------------------+-------------------------------------------------------------------------+
| -co/\\-\\-codon             | Codon codon-based trimming. *Default: off*                              |
+-----------------------------+-------------------------------------------------------------------------+
| -s/\\-\\-sequence           | Specifies sequence type of input file. *Default: auto-detect*           |
+-----------------------------+-------------------------------------------------------------------------+
| -if/\\-\\-input_file_format | Specify input file format*. *Default: auto-detect*                      |
+-----------------------------+-------------------------------------------------------------------------+
| -of/\\-\\-output_file_format| Specify output file format*. *Default: input file type*                 |
+-----------------------------+-------------------------------------------------------------------------+
| -l/\\-\\-log                | Create a log file. *Default: off*                                       |
+-----------------------------+-------------------------------------------------------------------------+
| -c/\\-\\-complementary      | Create a complementary alignment file. *Default: off*                   |
+-----------------------------+-------------------------------------------------------------------------+
| -a/\\-\\-auxiliary_file     | Auxiliary file. Currently used for specifying sites to trim in cst mode |
+-----------------------------+-------------------------------------------------------------------------+


\*Acceptable file formats include: 
`fasta <https://en.wikipedia.org/wiki/FASTA_format>`_,
`clustal <http://meme-suite.org/doc/clustalw-format.html>`_,
`maf <http://www.bx.psu.edu/~dcking/man/maf.xhtml>`_,
`mauve <http://darlinglab.org/mauve/user-guide/files.html>`_,
`phylip <http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html>`_,
`phylip-sequential <http://rosalind.info/glossary/phylip-format/>`_,
`phylip-relaxed <https://www.hiv.lanl.gov/content/sequence/FORMAT_CONVERSION/FormatExplain.html>`_,
`stockholm <https://en.wikipedia.org/wiki/Stockholm_format>`_


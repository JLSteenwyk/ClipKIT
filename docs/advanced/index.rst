Advanced Usage
==============

^^^^^

This section describes the various features and options of ClipKIT.

- Modes_
- Log_
- Complementary_
- `Miscellaneous options`_

.. _Modes:

Modes
-----

ClipKIT can be run with three different modes (kpi, gappy, kpi-gappy), which are specified with the -m/--mode argument.
*Default: 'gappy'*

* kpi will trim all sites that are not parsimony informative
* gappy will remove all sites that are above a threshold of gappyness (default: 0.9)
* kpi-gappy is the combination of kpi- and gappy-based trimming

.. code-block:: shell

	# kpi-based trimming
	clipkit <input> -m kpi

	# gappy-based trimming
	clipkit <input> -m gappy

	# kpi- and gappy-based trimming
	clipkit <input> -m kpi-gappy

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

.. _Complementary:

Complementary
-------------

Having an alignment of the sequences that were trimmed can be useful for other analyses. 
To obtain an alignment of the sequences that were trimmed, use the -c/\\-\\-complementary 
option.

.. code-block:: shell

	clipkit <input> -c

Output file with the suffix '.clipkit.complementary'

.. _`Miscellaneous options`:

Miscellaneous options
---------------------


+-----------------------------+-------------------------------------------------------------------+
| Option                      | Usage and meaning                                                 |
+=============================+===================================================================+
| -h/\\-\\-help               | Print help message                                                |
+-----------------------------+-------------------------------------------------------------------+
| -v/\\-\\-version            | Print software version                                            |
+-----------------------------+-------------------------------------------------------------------+
| -g/\\-\\-gaps     	      | Specify gappyness threshold (between 0 and 1). *Default: 0.9*     |
+-----------------------------+-------------------------------------------------------------------+
| -if/\\-\\-input_file_format | Specify input file format*. *Default: auto-detect*                |
+-----------------------------+-------------------------------------------------------------------+
| -of/\\-\\-input_file_format | Specify output file format*. *Default: input file type*           |
+-----------------------------+-------------------------------------------------------------------+

\*Acceptable file formats include: 
`fasta <https://en.wikipedia.org/wiki/FASTA_format>`_,
`clustal <http://meme-suite.org/doc/clustalw-format.html>`_,
`maf <http://www.bx.psu.edu/~dcking/man/maf.xhtml>`_,
`mauve <http://darlinglab.org/mauve/user-guide/files.html>`_,
`phylip <http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html>`_,
`phylip-sequential <http://rosalind.info/glossary/phylip-format/>`_,
`phylip-relaxed <https://www.hiv.lanl.gov/content/sequence/FORMAT_CONVERSION/FormatExplain.html>`_,
`stockholm <https://en.wikipedia.org/wiki/Stockholm_format>`_

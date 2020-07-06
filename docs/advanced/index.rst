.. image:: ../_static/img/ClipKIT_logo_top_only_v1.jpg
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/ClipKIT


Advanced Usage
==============

^^^^^

This section describes the various features and options of ClipKIT.

- Modes_
- Output_
- Log_
- Complementary_
- `All options`_

.. _Modes:

Modes
-----

ClipKIT can be run with five different modes (gappy, kpic, kpic-gappy, kpi, and kpi-gappy), which are specified with the -m/--mode argument.
*Default: 'gappy'*

* gappy: trim all sites that are above a threshold of gappyness (default: 0.9)
* kpic (alias: medium): keep only parismony informative and constant sites
* kpic-gappy (alias: medium-gappy): a combination of kpic- and gappy-based trimming
* kpi (alias: heavy): keep only parsimony informative sites
* kpi-gappy (alias: heavy-gappy): a combination of kpi- and gappy-based trimming

.. code-block:: shell

	# gappy-based trimming
	clipkit <input>
	clipkit <input> -m gappy

	# kpic-based trimming
	clipkit <input> -m kpic
	clipkit <input> -m medium

	# kpic- and gappy-based trimming
	clipkit <input> -m kpic-gappy
	clipkit <input> -m medium-gappy

	# kpi-based trimming
	clipkit <input> -m kpi
	clipkit <input> -m heavy

	# kpi- and gappy-based trimming
	clipkit <input> -m kpi-gappy
	clipkit <input> -m heavy-gappy

.. _Output:

Output
------

By default, output files will have the same name as the input file with the suffix ".clipkit"
appended to the name. Users can specify output file names with the -o option. 

.. code-block:: shell

	# specify output
	clipkit <input> -o <output>

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

.. _`All options`:

All options
---------------------


+-----------------------------+-------------------------------------------------------------------+
| Option                      | Usage and meaning                                                 |
+=============================+===================================================================+
| -h/\\-\\-help               | Print help message                                                |
+-----------------------------+-------------------------------------------------------------------+
| -v/\\-\\-version            | Print software version                                            |
+-----------------------------+-------------------------------------------------------------------+
| -o/\\-\\-output             | Specify output file name                                          |
+-----------------------------+-------------------------------------------------------------------+
| -g/\\-\\-gaps               | Specify gappyness threshold (between 0 and 1). *Default: 0.9*     |
+-----------------------------+-------------------------------------------------------------------+
| -if/\\-\\-input_file_format | Specify input file format*. *Default: auto-detect*                |
+-----------------------------+-------------------------------------------------------------------+
| -of/\\-\\-output_file_format| Specify output file format*. *Default: input file type*           |
+-----------------------------+-------------------------------------------------------------------+
| -l/\\-\\-log                | Create a log file. *Default: off*                                 |
+-----------------------------+-------------------------------------------------------------------+
| -c/--complementary          | Create a complementary alignment file. *Default: off*             |
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

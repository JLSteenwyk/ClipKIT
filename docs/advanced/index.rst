Advanced Usage
==============

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
- `Ends only`_
- Threads_
- `All options`_

|

.. _Modes:

Modes
-----

This section describes the trimming modes implemented in ClipKIT. If you are unsure which is appropriate for you,
**we recommend using the default smart-gap trimming mode**. 

ClipKIT can be run with eleven different modes, which are specified with the -m/--mode argument.
*Default: 'smart-gap'*

* smart-gap: dynamic determination of gaps threshold
* entropy: trim sites above a normalized Shannon entropy threshold (default: 0.8)
* gappy: trim all sites that are above a threshold of gappyness (default: 0.9)
* kpic: keep only parsimony informative and constant sites
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

	# entropy-based trimming
	clipkit <input> -m entropy

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
It can be useful to have information about each position in an alignment. For
example, this information could be used in alignment diagnostics, fine-tuning of trimming
parameters, etc. To create the log file, use the -l/\\-\\-log option. Using this option
will create a four-column file with the suffix 'clipkit.log'. *Default: off*

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

Output file with the suffix '.clipkit.complement'

|

.. _Codon:

Codon
-----

Trims codon-based alignments. If one position in a codon should be trimmed, the whole
codon will be trimmed. To conduct codon-based trimming, use the -co/\\-\\-codon argument.

.. code-block:: shell

	clipkit <input> --codon

    # or

	clipkit <input> -co

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
when using the kpi and kpic modes of trimming as well as an 
iteration of trimming that uses smart-gap. In entropy mode, this value
is treated as a normalized Shannon entropy threshold (default: 0.8).

To specify a gaps threshold, use the -g/\\-\\-gaps argument.

.. code-block:: shell

	clipkit <input> --gaps 0.4

    # or

	clipkit <input> -g 0.4

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

Use this option to specify that input sequences are amino acids.

.. code-block:: shell

	clipkit <input> -s nt

Use this option to specify that input sequences are nucleotides.

|

.. _`Ends only`:

Ends only
---------

For a given trimming mode, this option trims only sites at the ends of an alignment.
For example, if the sites that should be trimmed include
[0, 1, 2, 4, 5, 6, 14, 15, 16] for smart-gap mode and an alignment of
length 16, adding the ends_only mode will result in [0, 1, 2, 14, 15, 16]
being the sites that will be trimmed. Use this argument with -eo, \-\-ends_only.


.. code-block:: shell

	clipkit <input> -eo

	# or

	clipkit <input> --ends_only

|

.. _Threads:

Threads
-------

ClipKIT supports parallel processing for site classification and character frequency
calculations. For larger alignments, this can significantly speed up processing.

The number of threads can be specified using the -t/\\-\\-threads argument.
*Default: 1*

.. code-block:: shell

	# Single-threaded processing (default)
	clipkit <input>

	# Multi-threaded processing with 4 threads
	clipkit <input> -t 4

	# or
	clipkit <input> --threads 4

**Performance Notes:**

* Parallel processing is activated adaptively based on alignment size and requested threads
* For smaller alignments, single-threaded mode is typically faster due to multiprocessing overhead
* The optimal number of threads depends on your system and alignment size
* For KPI/KPIC family modes (kpi, kpi-gappy, kpi-smart-gap, kpic, kpic-gappy, kpic-smart-gap), ClipKIT may automatically use fewer threads than requested when that is expected to be faster
* Results are identical regardless of the number of threads used (fully reproducible)

|

.. _`Dry run`:

Dry run
-------

Use dry run mode to execute trimming and compute summary statistics without
writing alignment, complementary, or log output files.

.. code-block:: shell

	clipkit <input> --dry_run

|

.. _`Validate only`:

Validate only
-------------

Use validate-only mode to check input format and argument consistency (including
auxiliary file checks for ``cst`` mode) and then exit without trimming.

.. code-block:: shell

	clipkit <input> --validate_only

|

.. _`Report JSON`:

Report JSON
-----------

Write a machine-readable JSON report with run configuration and outcome details.

.. code-block:: shell

	# explicit report path
	clipkit <input> --report_json run_report.json

	# default report path: <output>.report.json
	clipkit <input> --report_json

|

.. _`All options`:

All options
---------------------


.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Option
     - Usage and meaning
   * - ``-h/--help``
     - Print help message.
   * - ``-v/--version``
     - Print software version.
   * - ``-m/--mode``
     - Specify trimming mode (including ``entropy``). *Default: smart-gap*.
   * - ``-o/--output``
     - Specify output file name.
   * - ``-g/--gaps``
     - Specify gappyness threshold (between 0 and 1). *Default: 0.9*.
   * - ``-gc/--gap_characters``
     - Specify gap characters used in input file (AAs: ``Xx-?*``; NTs: ``XxNn-?*``).
   * - ``-co/--codon``
     - Conduct codon-based trimming. *Default: off*.
   * - ``-s/--sequence_type``
     - Specify sequence type of input file (``aa`` or ``nt``). *Default: auto-detect*.
   * - ``-if/--input_file_format``
     - Specify input file format*. *Default: auto-detect*.
   * - ``-of/--output_file_format``
     - Specify output file format*. *Default: input file type*.
   * - ``-l/--log``
     - Create a log file. *Default: off*.
   * - ``-c/--complementary``
     - Create a complementary alignment file. *Default: off*.
   * - ``-a/--auxiliary_file``
     - Auxiliary file used for specifying sites to trim in ``cst`` mode.
   * - ``-eo/--ends_only``
     - Trim only sites at alignment ends that would otherwise be removed. *Default: off*.
   * - ``-q/--quiet``
     - Disable logging to stdout. *Default: off*.
   * - ``-t/--threads``
     - Requested threads for parallel processing; KPI/KPIC modes may auto-tune lower. *Default: 1*.
   * - ``--dry_run``
     - Run trimming/stat calculations but skip writing output files. *Default: off*.
   * - ``--validate_only``
     - Validate inputs/arguments and exit without trimming. *Default: off*.
   * - ``--report_json [path]``
     - Write a JSON run report; if no path is given, uses ``<output>.report.json``.

\*Acceptable file formats include: 
`fasta <https://en.wikipedia.org/wiki/FASTA_format>`_,
`clustal <http://meme-suite.org/doc/clustalw-format.html>`_,
`maf <http://www.bx.psu.edu/~dcking/man/maf.xhtml>`_,
`mauve <http://darlinglab.org/mauve/user-guide/files.html>`_,
`phylip <http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html>`_,
`phylip-sequential <http://rosalind.info/glossary/phylip-format/>`_,
`phylip-relaxed <https://www.hiv.lanl.gov/content/sequence/FORMAT_CONVERSION/FormatExplain.html>`_,
`stockholm <https://en.wikipedia.org/wiki/Stockholm_format>`_

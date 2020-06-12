.. image:: _static/img/ClipKIT_logo_top_only_v1.jpg
   :width: 55%
   :align: center

^^^^^


ClipKIT is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes those that display characteristics poor phylogenetic signal.


If you found clipkit useful, please cite *ClipKIT: a multiple sequence alignment-trimming algorithm for accurate phylogenomic inference*. bioRxiv. doi: |doiLink|_.

.. _doiLink: https://www.biorxiv.org/content/10.1101/2020.06.08.140384v1
.. |doiLink| replace:: 10.1101/2020.06.08.140384 


Quick Start
-----------
**1) Installation**

To install, use the following commands:

.. code-block:: shell

	pip install clipkit

|

To install from source, use the following commands:

.. code-block:: shell

	git clone https://github.com/JLSteenwyk/ClipKIT.git
	cd ClipKIT/
	make install

If you run into permission errors when executing *make install*, create a 
virtual environemnt for your installation:

.. code-block:: shell

	git clone https://github.com/JLSteenwyk/ClipKIT.git
	cd ClipKIT/
	python -m venv .venv
	source .venv/bin/activate
	make install

Note, the virtual environment must be activated to use clipkit.

|

**2) Usage**

To use ClipKIT in its simpliest form, execute the following command:

.. code-block:: shell

	clipkit <input>

Output file with the suffix ".clipkit"

|

^^^^

.. toctree::
	:maxdepth: 4

	about/index
	advanced/index
	performance_assessment/index

^^^^


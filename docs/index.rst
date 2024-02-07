.. image:: _static/img/ClipKIT_logo_top_only_v1.jpg
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/ClipKIT

^^^^^


ClipKIT is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes those that display characteristics poor phylogenetic signal.


If you found clipkit useful, please cite *ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference*. Steenwyk et al. 2020, PLoS Biology. doi: |doiLink|_.

.. _doiLink: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007
.. |doiLink| replace:: 10.1371/journal.pbio.3001007


Quick Start
-----------
These two lines represent the simpliest method to rapidly install and run ClipKIT.

.. code-block:: shell

	# install
	pip install clipkit
	# run
	clipkit input.fa

Below are more detailed instructions, including alternative installation methods.

**1) Installation**

To help ensure ClipKIT can be installed using your favorite workflow, we have made ClipKIT available from pip, source, and the anaconda cloud.

**Install from pip**

To install from pip, use the following commands:

.. code-block:: shell

	# create virtual environment
	python -m venv venv
	# activate virtual environment
	source venv/bin/activate
	# install clipkit
	pip install clipkit

**Note, the virtual environment must be activated to use clipkit.**

|

**Install from source**

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the following commands:

.. code-block:: shell

	# download
	git clone https://github.com/JLSteenwyk/ClipKIT.git
	cd ClipKIT/
	# create virtual environment
	python -m venv venv
	# activate virtual environment
	source venv/bin/activate
	# install
	make install

To deactivate your virtual environment, use the following command:

.. code-block:: shell

	# deactivate virtual environment
	deactivate

**Note, the virtual environment must be activated to use clipkit.**

.. code-block:: shell

	git clone https://github.com/JLSteenwyk/ClipKIT.git
	cd ClipKIT/
	make install

If you run into permission errors when executing *make install*, create a 
virtual environemnt for your installation:

.. code-block:: shell

	git clone https://github.com/JLSteenwyk/ClipKIT.git
	cd ClipKIT/
	python -m venv venv 
	source venv/bin/activate
	make install

Note, the virtual environment must be activated to use clipkit.

|

**Install from anaconda**

To install via anaconda, execute the following command:

.. code-block:: shell

	conda install bioconda::clipkit

Visit here for more information: https://anaconda.org/bioconda/clipkit

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
	change_log/index
	other_software/index
	frequently_asked_questions/index

^^^^


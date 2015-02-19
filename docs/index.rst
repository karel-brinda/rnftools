.. MISmash documentation master file, created by
   sphinx-quickstart on Mon Feb  9 11:46:57 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MIShmash's documentation!
====================================

MIShmash is a pipeline for simulation of Next-Generation Sequencing reads. It simulates reads using
existing read simulators and combines the obtained reads into a single set. The pipeline is based
on SnakeMake. All required software is installed when requested.

Prerequisities
--------------

- Unix-like operating system (Linux, MacOS, etc.)
- Python 3.2 -- if not installed yet, Anaconda (http://continuum.io/downloads) is a recommended distribution
- SnakeMake (see http://bitbucket.org/johanneskoester/snakemake/), usually can be installed or upgraded using
	.. code-block:: bash
		pip install snakemake


Installation / upgrade to the latest version
--------------------------------------------

- using PIP:

	.. code-block:: bash
	
		pip install --upgrade mishmash

- from GIT (the last development version):

	.. code-block:: bash
	
		git clone http://github.com/karel-brinda/mishmash
		cd mishmash
		python3 setup.py install


Usage
-----

Create a directory where you want to simulate reads.


Supported read simulators
-------------------------

ART Illumina
~~~~~~~~~~~~

Syntax: 

.. code-block:: python

	mishmash.ArtIllumina(
		fa
		coverate=0,
		number_of_reads=0,
		read_length_1=100,
		read_length_2=0,
		other_params="",
		distance=500,
		distance_deviation=50.0,
		rng_seed=1,
	)


DwgSim
~~~~~~

Syntax:

.. code-block:: python

	mishmash.DwgSim(
		fa,
		coverage=0,
		number_of_reads=0,
		read_length_1=100,
		read_length_2=0,
		other_params="",
		distance=500,
		distance_deviation=50.0,
		rng_seed=1,
	)


Contents:

.. toctree::
   :maxdepth: 2

.. automodule:: mishmash
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


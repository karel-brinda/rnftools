.. MISmash documentation master file, created by
   sphinx-quickstart on Mon Feb  9 11:46:57 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================================
Welcome to MIShmash's documentation!
====================================

MIShmash is a pipeline for simulation of Next-Generation Sequencing reads. It simulates reads using
existing read simulators and combines the obtained reads into a single set. The pipeline is based
on SnakeMake. All required software is installed when requested.

.. contents:: Table of Contents
   :depth: 2


Prerequisities
==============

- **Unix-like operating system** (Linux, MacOS, etc.)
- **Python 3.2 or higher**
	- if not installed yet, Anaconda (http://continuum.io/downloads) is a recommended distribution
- **PIP**
	- for installation, see https://pip.pypa.io/en/latest/installing.html
- **SnakeMake** - see http://bitbucket.org/johanneskoester/snakemake/)
	- it can be usually installed using

		.. code-block:: bash
		
			pip install snakemake


Installation / upgrade to the latest version
============================================

- using PIP (recommended):

	.. code-block:: bash
	
		pip install --upgrade mishmash

- from GIT (the last development version):

	.. code-block:: bash
	
		git clone http://github.com/karel-brinda/mishmash
		cd mishmash
		python3 setup.py install


Usage
=====

Create a directory where you want to simulate reads.


Supported read simulators
=========================

Explanation of shared parameters:

- ``fa`` -- reference (FASTA file)
- ``coverage`` -- average coverage (0 = non-specified)
- ``number_of_reads`` -- number of reads (0 = non-specified)
- ``read_length_1`` -- length of the first end of a read
- ``read_length_2`` -- lenghh of the second end of a read (0 => single-end simulation)
- ``other_params`` -- other parameters for the given simulator (shell string)
- ``distance`` -- mean inner distance between ends of a read
- ``distance_deviation`` -- its deviation
- ``rng_seed`` -- seed for random number generator

Remarks:

- ``coverage`` or ``number_of_reads`` must be equal to zero


ART Illumina
------------

Example:
~~~~~~~~

.. code-block:: python

	mishmash.ArtIllumina(
		fa="my_reference.fa",
		number_of_reads=10000,
		read_length_1=100,
		read_length_2=0,
	)


Syntax: 
~~~~~~~

.. code-block:: python

	mishmash.ArtIllumina(
		fa
		coverage=0,
		number_of_reads=0,
		read_length_1=100,
		read_length_2=0,
		other_params="",
		distance=500,
		distance_deviation=50.0,
		rng_seed=1,
	)


DwgSim
------

Example:
~~~~~~~~

.. code-block:: python

	mishmash.DwgSim(
		fa="my_referenc.fa",
		number_of_reads=10000,
		read_length_1=100,
		read_length_2=100,
	)

Syntax:
~~~~~~~

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

Remarks:
~~~~~~~~

- for pair-end read simulation, ``read_length_1`` must equal to ``read_length_2``



Contents
========

.. toctree::
   :maxdepth: 2

.. automodule:: mishmash
   :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


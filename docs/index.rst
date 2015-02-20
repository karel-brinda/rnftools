.. MIShmash documentation master file, created by
   sphinx-quickstart on Mon Feb  9 11:46:57 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MIShmash
========

MIShmash is a pipeline for simulation of Next-Generation Sequencing reads. It simulates reads using
existing read simulators and combines the obtained reads into a single set. The pipeline is based
on SnakeMake. All required software is installed automatically when requested.

Table of Contents
-----------------

.. contents::
   :depth: 2


Requirements
------------

* Unix-like operating system (Linux, MacOS, etc.).
* `Python`_ 3.2+. If not installed yet, Anaconda is a recommended Python distribution.
* `SnakeMake`_
* `PIP`_

.. _Python: http://python.org
.. _Anaconda: http://continuum.io/downloads
.. _SnakeMake: http://bitbucket.org/johanneskoester/snakemake/
.. _PIP: http://pip.pypa.io/en/latest/installing.html

Installation
------------

* using PIP (recommended):

	.. code-block:: bash
	
		pip install --upgrade mishmash

* from GIT (the last development version):

	.. code-block:: bash
	
		git clone http://github.com/karel-brinda/mishmash
		cd mishmash
		python3 setup.py install


Usage
-----

Technically, MIShmash works as a SnakeMake pipeline. To use it, create an empty directory where the reads will be simulated.
Create there an empty file named ``Snakefile``, which will serve as a configuration script.
Then save the following content into it:

.. code-block:: python
	
	# required line, it should be the first line of all your configuration scripts
	import mishmash

	# this line tells MIShmash that there will be a new sample
	mishmash.sample("new_sample")
	
	# then you can add arbitrary number of sources
	mishmash.ArtIllumina(
		fa="my_fast.fa",
		number_of_reads=10000,
		read_length_1=100,
		read_length_2=0,
	)
	
	# if you want to create more simulated samples, call again the mishmash.sample
	# function but with another sample name


	# these lines are mandatory as the last lines of the file
	include: mishmash.include
	rule: input: mishmash.output

Supported read simulators
-------------------------

Explanation of the shared parameters:

+----------------------------+--------------------------------------------------------------------+
| ``fa``                     | reference (FASTA file)                                             |
+----------------------------+--------------------------------------------------------------------+
| ``coverage``               | average coverage (0 = non-specified)                               |
+----------------------------+--------------------------------------------------------------------+
| ``number_of_reads``        | number of reads (0 = non-specified)                                |
+----------------------------+--------------------------------------------------------------------+
| ``read_length_1``          | length of the first end of a read                                  |
+----------------------------+--------------------------------------------------------------------+
| ``read_length_2``          | length of the second end of a read (0 => single-end simulation)    |
+----------------------------+--------------------------------------------------------------------+
| ``other_params``           | other parameters for the given simulator (shell string)            |
+----------------------------+--------------------------------------------------------------------+
| ``distance``               | mean inner distance between ends of a read                         |
+----------------------------+--------------------------------------------------------------------+
| ``distance_deviation``     | its deviation                                                      |
+----------------------------+--------------------------------------------------------------------+
| ``rng_seed``               | seed for random number generator                                   |
+----------------------------+--------------------------------------------------------------------+

Remarks:

* ``coverage`` or ``number_of_reads`` must be equal to zero


ART Illumina
^^^^^^^^^^^^

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
^^^^^^

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

* for pair-end read simulation, ``read_length_1`` must equal to ``read_length_2``



Contents
--------

.. toctree::
   :maxdepth: 2

.. automodule:: mishmash
   :members:


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


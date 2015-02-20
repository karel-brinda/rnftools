.. LAVEnder documentation master file, created by
   sphinx-quickstart on Fri Feb 20 11:33:55 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LAVEnder
========

LAVEnder is a RNF-compatible evaluation tool for Next-Generation Sequencing mappers.


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
	
		pip install --upgrade lavender

* from GIT (the last development version):

	.. code-block:: bash
	
		git clone http://github.com/karel-brinda/lavender
		cd lavender
		python3 setup.py install


Usage
-----

Technically, MIShmash works as a SnakeMake pipeline.

1. Create simulated reads using a RNF-compatible read simulator, for example `MIShmash`_
2. Map the simulated reads to reference and obtained BAM file put into a standalone directory.
3. Create an empty directory where a report will be created.
4. Create there an empty file named ``Snakefile``, which will serve as a configuration script. Then save the following content into it:

	.. code-block:: python

		# required line, it should be the first line of all you LAVEnder scripts
		import lavender

		#
		lavender.Report(
			bam_dirs=["my_bam_directory"],
			name="My amazing report"
		)

		# these lines are mandatory last lines of the file
		include:lavender.include
		rule: input: lavender.input

5. Run ``snakemake`` in the directory


.. _MIShmash: http://mishmash.rtfd.org



Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


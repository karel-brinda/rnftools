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
* `Python`_ 3.2+. If not installed yet, `Anaconda`_ is a recommended Python distribution.
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

Technically, MIShmash works as a `SnakeMake`_ pipeline.

1. Create simulated reads using a RNF-compatible read simulator, for example `MIShmash`_
2. Map the simulated reads to your reference with your tested mappers and put the obtained BAM files into a standalone directory.
3. Create an empty directory where a report will be created.
4. Create there a file named ``Snakefile``, which will serve as a configuration script. Save the following content into it and adjust it as you want.

	.. code-block:: python

		# a mandatory line
		import lavender

		# this command says: "create a report called 'my_amazing_report' for BAM files
		# in directory 'my_bam_directory'"
		lavender.Report(
			bam_dirs=["my_bam_directory"],
			name="My amazing report"
		)

		# two mandatory lines
		include: lavender.include
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

.. automodule:: lavender

.. autoclass:: lavender.Report
	:members:
	:inherited-members:
	:show-inheritance:

.. autoclass:: lavender.Pannel
	:members:
	:inherited-members:
	:show-inheritance:

.. autoclass:: lavender.Bam
	:members:
	:inherited-members:
	:show-inheritance:


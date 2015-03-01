Mapper evaluation
=================

.. contents::
   :depth: 2


LAVender
--------

LAVEnder is a RNF-compatible evaluation tool for evaluation of mappers of Next-Generation Sequencing reads.



Usage
^^^^^

1. Create simulated reads using a RNF-compatible read simulator, for example `MIShmash`_
2. Map the simulated reads to your reference with your tested mappers and put the obtained BAM files into a standalone directory.
3. Create an empty directory where a report will be created.
4. Create there a file named ``Snakefile``, which will serve as a configuration script. Save the following content into it and adjust it as you want.

	.. code-block:: python

		# a mandatory line
		import rnftools.lavender

		# this command says: "create a report called 'my_amazing_report' for BAM files
		# in directory 'my_bam_directory'"
		lavender.Report(
			bam_dirs=["my_bam_directory"],
			name="My amazing report"
		)

		# two mandatory lines
		include: lavender.include()
		rule: input: rnftools.input()

5. Run ``snakemake`` in the directory



Technical documentation
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: rnftools.lavender.Report
	:members:
	:inherited-members:
	:show-inheritance:

.. autoclass:: rnftools.lavender.Panel
	:members:
	:inherited-members:
	:show-inheritance:

.. autoclass:: rnftools.lavender.Bam
	:members:
	:inherited-members:
	:show-inheritance:


Simulation of NGS reads
======================

.. contents::
   :depth: 3


MIShmash
--------

MIShmash is a part of RNFtools responsible for simulation of Next-Generation Sequencing reads. It employs
existing read simulators and combines the obtained reads into bigger single set.


Usage
-----

Create there an empty file named ``Snakefile``, which will serve as a configuration script.
Then save the following content into it:

.. code-block:: python
	
	# required line, it should be the first line of all your configuration scripts
	import rnftools

	# this line tells MIShmash that there will be a new sample
	rnftools.mishmash.sample("new_sample")
	
	# then you can add arbitrary number of sources
	rnftools.mishmash.ArtIllumina(
		fa="my_fast.fa",
		number_of_reads=10000,
		read_length_1=100,
		read_length_2=0,
	)
	
	# if you want to create more simulated samples, call again the mishmash.sample
	# function but with another sample name


	# these lines are mandatory as the last lines of the file
	include: rnftools.mishmash.include()
	rule: input: rnftools.mishmash.output()


Supported simulators
--------------------

ART
^^^

+----------------+----------------------------------------------------------------------------------------------------------+
| URL:           | http://www.niehs.nih.gov/research/resources/software/biostatistics/art/                                  |
+----------------+----------------------------------------------------------------------------------------------------------+
| Publication:   | Huang, W. _et al._                                                                                       |
|                | ART: a next-generation sequencing read simulator.                                                        |
|                | _Bioinformatics_ *28*\(4), pp. 593--594, 2011.                                                           |
+----------------+----------------------------------------------------------------------------------------------------------+

ART Illumina
~~~~~~~~~~~~

.. autoclass:: rnftools.mishmash.ArtIllumina
        :members:
        :inherited-members:
        :show-inheritance:

CuReSim
^^^^^^^

.. autoclass:: rnftools.mishmash.CuReSim
        :members:
        :inherited-members:
        :show-inheritance:

DwgSim
^^^^^^

.. autoclass:: rnftools.mishmash.DwgSim
        :members:
        :inherited-members:
        :show-inheritance:

WgSim
^^^^^

.. autoclass:: rnftools.mishmash.WgSim
        :members:
        :inherited-members:
        :show-inheritance:




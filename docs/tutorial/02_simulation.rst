Read simulation
===============

In this chapter, we show on several basic example how to simulate reads using a
component of RNFtools called MISmash.


Basic example
-------------

First, let us show how to simulate reads from a single genome, stored in a FASTA
file, using a single simulator. A corresponding RNFtools configuration
script can look as follows:

.. literalinclude:: ../../examples/01_tutorial/02_simulation/01_simple_simulation/Snakefile
	:language: python
	:linenos:


Lines 1, 15, 16 have been already described in the previous chapter. Function
``rnf.mishmash.sample`` (line 3) initializes a new sample of simulated
"single-end reads". When a new instance of ``rnftools.mishmash.ArtIllumina`` is
created (line 8), it is automatically registered to the last created sample.
This class is used for simulating reads by the Art Illumina read simulator. The
parameters signalize parameters of the simulation: ``fasta`` is the reference
file, ``number_of_read_tuples`` sets number of simulated read tuples, and
``read_length_1`` and ``read_length_2`` indicate lengths of simulated reads.

Within the RNF framework, a single simulated unit is a read tuple, which
consists of one or more reads. For more details, see the `RNFtools paper`_.
``read_length_2=0`` implies "single-end" simulation (in our notation: a single
read in every read tuple).

When we run ``snakemake``, reads are simulated and we obtain the final
``simple_example.fq`` file with all the simulated reads.

RNFtools supports several different read simulators. Their use is similar,
though their interfaces are slightly different. A full documentation of all the
supported simulators with all their parameters is available on the page
:ref:`mishmash`.


Simulation of 'paired-end' reads
--------------------------------

To simulate "paired-end reads" (i.e., read tuples of two reads), two minor
changes must be made in the original ``Snakefile``. First,
``rnftools.mishmash.sample`` must be called with ``reads_in_tuple=2``, then
the length of second reads of a tuple must be set to a non-zero value.

Then the final ``Snakefile`` can be:

.. literalinclude:: ../../examples/01_tutorial/02_simulation/02_paired_end_reads/Snakefile
	:language: python
	:linenos:


Different simulator
-------------------

To change a simulator in our example, we just replace
``rnftools.mishmash.ArtIllumina`` by another simulator, e.g.,
``rnftools.mishmash.ArtIllumina``. Parameters as ``fasta``, ``read_length_1``,
``read_length_2``, and ``number_of_read_tuples`` are the same for all the
simulators, but with several limitations:

* CuReSim supports only "single-end reads".
* ART Illumina in "paired-end" mode can simulate only reads of equal lengths.


.. literalinclude:: ../../examples/01_tutorial/02_simulation/03_different_simulator/Snakefile
	:language: python
	:linenos:


More genomes
------------

To simulate reads from multiples reference within a single sample (in order to
simulate, e.g., a metagenome or a contamination), we create a new instance of
class of a simulator for each reference.

In the example below, we are simulating 10.000 read tuples from two
reference genomes.

.. literalinclude:: ../../examples/01_tutorial/02_simulation/04_more_genomes/Snakefile
	:language: python
	:linenos:

Once reads are simulated for each of the references, they are mixed and put
into the resulting FASTQ file.


More samples
------------

We can also simulated multiple samples using a single ``Snakemake``.

.. literalinclude:: ../../examples/01_tutorial/02_simulation/05_more_samples/Snakefile
	:language: python
	:linenos:



Sequence extraction
-------------------

It may be sometimes useful to extract certain sequences from the reference file
before the simulation itself. For instance, reads from each sequence could be
simulated with a different coverage. For this purpose, we can use the
``sequences`` parameter. Sequences for extraction can be specified either by
their number in the original FASTA file (starting from 0), or using their name.

.. literalinclude:: ../../examples/01_tutorial/02_simulation/06_sequence_extraction/Snakefile
	:language: python
	:linenos:

Non-standard parameters
-----------------------

Not all command-line parameters of every simulator are directly supported by
RNFtools. However, such parameters can still be passed through the parameter
``other_params`` like in this example:

.. literalinclude:: ../../examples/01_tutorial/02_simulation/07_nonstandard_parameters/Snakefile
	:language: python
	:linenos:

.. _`RNFtools paper`: http://dx.doi.org/10.1093/bioinformatics/btv524

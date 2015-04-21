Read simulation
---------------

In this chapter we show how to simulate reads in basic simulation. For this task, component called MIShmash is used.


Basic example
"""""""""""""

First let us show how to simulate reads from a single genome (saved in a FASTA file) using a single simulator. Then a
corresponding RNFtools configuration script can look as follows:

.. literalinclude:: ../../examples/tutorial/02_simulation/01_simple_simulation/Snakefile
	:language: python
	:linenos:


Lines 1, 2, 13, 14 were already described in the previous chapter. Function ``rnf.mishmash.sample`` on line 4 initializes a new sample of simulated "single-end reads". When a new instance of ``rnftools.mishmash.ArtIllumina`` is created on line 6, it is automatically registered to the last created sample. This class is used for simulating reads by Art Illumina read simulator. The parameters signalize parameters of the simulation: ``fasta`` is the reference file, ``number_of_read_tuples`` sets number of simulated read tuples, and  ``read_length_1`` and ``read_length_2`` indicate lengths of simulated reads.

Within RNF framework, a single simulated unit is a read tuple which consists of one or more reads. For more details, see the RNFtools paper. ``read_length_2=0`` implies "single-end" simulation (in our notation: single read in every read tuple).

In our code, we use also SMBL library (import on line 2, usage on line 7) for automatic download of an example of a reference genome (SMBL library is also internally used by RNFtools for automatic installation of all employed programs). You can change the ``fasta`` parameter to some existing FASTA file, e.g., ``fasta="reference.fa"``.

Let us run ``snakemake``. All employed programs are automatically downloaded and compiled (by the SMBL library).
An example of a reference is downloaded. Finally reads are simulated and you obtain final ``simple_example.fq`` file with simulated reads.

All programs were installed into ``~/.smbl/bin/`` and the example of a FASTA file to ``~/.smbl/fa/``. These are default output directories for the SMBL library.

RNFtools can work with more read simulators and usage of all of them is very similar. However there exist some differences between their interfaces. Full documentation of simulators with all parameters is available in Reference/MIShmash.


Simulation of 'paired-end' reads
""""""""""""""""""""""""""""""""

To simulate "paired-end reads" (i.e., read tuples of two reads), two minor changes must be done in the original ``Snakefile``. First, ``rnftools.mishmash.sample`` must be called with ``reads_in_tuple=2``, then length of second read of a tuple must be set to a non-zero value.

Then the final ``Snakefile`` can be:

.. literalinclude:: ../../examples/tutorial/02_simulation/02_paired_end_reads/Snakefile
	:language: python
	:linenos:


Different simulator
"""""""""""""""""""

To change simulator in our example, just replace ``rnftools.mishmash.ArtIllumina`` by class of another simulator, e.g., ``rnftools.mishmash.ArtIllumina``. Parameters as ``fasta``, ``read_length_1``, ``read_length_2``, and ``number_of_read_tuples`` are shared by all of them.

When you are changing the used simulator, be aware of these limitations:

* CuReSim supports only "single-end reads".
* ART Illumina in "paired-end" mode can simulate only reads of equal lengths.


.. literalinclude:: ../../examples/tutorial/02_simulation/03_different_simulator/Snakefile
	:language: python
	:linenos:


More genomes
""""""""""""

To simulate reads from more genomes and mix them in one sample (in order to simulate, e.g., metagenome or contamination), create a new instance of class of a simulator. In this example, we are simulating reads from two examples of reference genomes
and from both of them, we simulate 10.000 read tuples. 

.. literalinclude:: ../../examples/tutorial/02_simulation/04_more_genomes/Snakefile
	:language: python
	:linenos:


More samples
""""""""""""

To create more samples, one has just to call ``rnftools.mishmash.sample`` more times in the file. In the following example, we create two simulations, one with "single-end reads", the other one with "paired-end reads".

.. literalinclude:: ../../examples/tutorial/02_simulation/05_more_samples/Snakefile
	:language: python
	:linenos:


Nonstandard parameters
""""""""""""""""""""""

Some specific parameter of a read simulator may not be supported by the corresponding class. Such parameters can be used anyway since there, for all simulators' classes, exists a parameter ``other_params``. 

.. literalinclude:: ../../examples/tutorial/02_simulation/06_nonstandard_parameters/Snakefile
	:language: python
	:linenos:

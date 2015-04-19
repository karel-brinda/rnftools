Read simulation
---------------

In this chapter we show how to simulate reads in basic simulation. For this task, component called MIShmash is used.


Simple simulation
"""""""""""""""""

First let us show how to simulate reads from a single genome (saved in a reference file) using a single simulator. Then the
corresponding RNFtools configuration script looks as follows:

.. literalinclude:: ../../examples/tutorial/02_simulation/01_simple_simulation/Snakefile
	:language: python
	:linenos:


Line starting 4 initializes a batch of simulated reads. Creating the
``rnftools.mishmash.ArtIllumina`` object signalizes that reads will be simualted by ``art_illumina`` from reference ``fasta`` with given number of read tuples and of given read lengths (``read_length_1`` and ``read_length_2``). In our notation, a single simulated unit is a read tuple which consists of`one or more reads. For more details, see the RNFtools paper. ``read_length_2=0`` implies "single-end" simulation (in our notation: single read in every read tuple).

In our code, we use also SMBL library to download automatically a reference (this library is also internally used by RNFtools). You can change the ``fasta`` parameter to some existing FASTA file, e.g., ``fasta="reference.fa"``.

When you run ``Snakemake``, more things happen. All employed programs are automatically downloaded and compiled (by the SMBL library).
An example reference is downloaded. Finally reads are simulated and you obtain final ``simple_example.fq`` file with simulated reads.

All programs were installed into ``~/.smbl/bin/`` and the example FASTA file to ``~/.smbl/fa/``.

RNFtools can work with more read simulators and work with them is very similar, however there exist differences between
their interfaces. Full documentation of simulators with all parameters is available in Reference/MIShmash.


Simulation of 'paired-end' reads
""""""""""""""""""""""""""""""""

.. literalinclude:: ../../examples/tutorial/02_simulation/02_paired_end_reads/Snakefile
	:language: python
	:linenos:


Different simulator
"""""""""""""""""""

To change simulator in our example, just replace ``rnftools.mishmash.ArtIllumina`` by class of another simulator, e.g., ```rnftools.mishmash.ArtIllumina```. Parameters like ``fasta``, ``read_length_1``, ``read_length_2``, and ``number_of_read_tuples`` are same for all of them.

When you are changing the employed simulator, be aware of these limitations:

* CuReSim supports only single-end reads.
* ART Illumina in paired-end mode can simulate only reads of equal lengths.

.. literalinclude:: ../../examples/tutorial/02_simulation/03_different_simulator/Snakefile
	:language: python
	:linenos:


More genomes
""""""""""""

To simulate reads from more genomes and mix them in one sample (in order to simulate, e.g., metagenome or contamination), add a new 

.. literalinclude:: ../../examples/tutorial/02_simulation/04_more_genomes/Snakefile
	:language: python
	:linenos:


More samples
""""""""""""

.. literalinclude:: ../../examples/tutorial/02_simulation/05_more_samples/Snakefile
	:language: python
	:linenos:



Nonstandard parameters
""""""""""""""""""""""

.. literalinclude:: ../../examples/tutorial/02_simulation/06_nonstandard_parameters/Snakefile
	:language: python
	:linenos:


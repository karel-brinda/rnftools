Simulation
----------

In this chapter we show how to simulate reads in basic simulation. For this task, component called MIShmash is used.


Simple simulation
"""""""""""""""""

First let us show how to simulate reads from a single genome (saved in a reference file) using a single simulator. Then the
corresponding RNFtools configuration script looks as follows:

.. code-block:: python

	import rnftools
	import smbl

	rnftools.mishmash.sample("simple_example",reads_in_tuple=1)

	rnftools.mishmash.ArtIllumina(
		fasta=smbl.fasta.EXAMPLE,
		number_of_read_tuples=10000,
		read_length_1=100,
		read_length_2=0,
	)

	include: rnftools.include()
	rule: input: rnftools.input()


Line starting with ``rnftools.mishmash.sample`` initializes a batch of simulated reads. Creating the
``rnftools.mishmash.ArtIllumina`` object signalizes that reads will be simualted by ``art_illumina`` from reference ``fasta`` with given number of read tuples and of given read lengths (``read_length_1`` and ``read_length_2``). In our notation, a single simulated unit is a read tuple which consists of`one or more reads. For more details, see the RNFtools paper. ``read_length_2=0`` implies "single-end" simulation (in our notation: single read in every read tuple).

In our code, we use also SMBL library to download automatically a reference (this library is also internally used by RNFtools). You can change the ``fasta`` parameter to some existing FASTA file, e.g., ``fasta="reference.fa"``.

When you run ``Snakemake``, more things happen. All employed programs are automatically downloaded and compiled (by the SMBL library).
An example reference is downloaded. Finally reads are simulated and you obtain final ``simple_example.fq`` file with simulated reads.

All programs were installed into ``~/.smbl/bin/`` and the example FASTA file to ``~/.smbl/fa/``.

Simulator change
""""""""""""""""

More genomes
""""""""""""

More samples
""""""""""""




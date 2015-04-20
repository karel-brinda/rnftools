Hello world!
------------

In this chapter we show how to create and use RNFtools on a Hello world example. Little knowledge
of Python can be helpful but it is not required.

RNFtools is based on Snakemake, a Python-based Make-like build system. You create
simple configuration Python scripts and RNFtools subsequently creates a set of rules to be run by Snakemake.
The rules can be then executed in a single thread, in parallel (``--cores`` Snakemake parameter), or on a
cluster. For more details about Snakemake, please see its documentation.

This approach enables creating big and reproducible pipelines, easy for sharing (you just publish
you configuration script).

Every RNFtools script consists of three parts:

.. literalinclude:: ../../examples/tutorial/01_first_snakefile/Snakefile
	:language: python
	:linenos:


As it has been mentioned in comment, all your code must be inserted into part 2. Now save this file with name ``Snakefile`` and run

.. code-block:: bash

	snakemake

in that directory. Nothing will happed, only your "Hello world" message and few informative messages will
be printed.

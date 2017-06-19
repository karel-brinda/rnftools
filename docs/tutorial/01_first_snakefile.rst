Hello world!
------------

In this chapter, we show how to create and use RNFtools on a Hello world
example. Little knowledge of Python can be helpful, but it is not required.

RNFtools is based on Snakemake_, a Python-based Make-like build system. To
simulate reads or evaluate alignments, you create simple configuration Python
scripts and RNFtools subsequently creates a set of rules to be run by
Snakemake.  The rules can be then executed in a single thread, in parallel
(``--cores`` Snakemake parameter), or on a cluster.  For more details about
Snakemake, please see its `documentation <https://snakemake.readthedocs.io>`_.

This approach allows to create big and reproducible pipelines, which are easy
to share (it suffices to publish a single configuration script).

Every RNFtools script consists of three parts:

.. literalinclude:: ../../examples/01_tutorial/01_first_snakefile/Snakefile
	:language: python
	:linenos:


As it is mentioned in the second comment, all your code should be inserted into
part 2. Now save this file with name ``Snakefile`` and run

.. code-block:: bash

	snakemake

in the same directory. A "Hello world" message will be printed, together with
several Snakemake informative messages.

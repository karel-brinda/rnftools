First Snakefile
---------------

RNFtools is in fact a Snakemake package. Therefore, configuration scripts
use Python language inside. For every task, a RNFtools script must be created.

Every RNFtools script consists of three parts:

.. code-block:: python

	# 1) a mandatory line loading rnftools package
	import rnftools

	# 2) your code
	print("Hello world")

	# 3) two lines saying that RNFtools rules files must be included and setting an initial rule
	include: rnftools.include()
	rule: input: rnftools.input()


All our code must be inserted into part 2. Now save this file with name ``Snakefile`` and run

.. code-block:: bash

	snakemake

in that directory. Nothing will happed, only your "Hello world" and few informative messages will
be printed.


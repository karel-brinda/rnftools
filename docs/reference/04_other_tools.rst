Other tools
===========

``rnf-join-fq.py``
------------------

A program for joining FASTQ files with reads in RNF format.

.. code::

	$ rnf-join-fq.py  -h
	usage: rnf-join-fq.py [-h] -i inp [inp ...] -m mode -o out

	Join FASTQ files with reads in RNF format.

	optional arguments:
	  -h, --help        show this help message and exit
	  -i inp [inp ...]  input FASTQ files
	  -m mode           mode for joining files (single-end / paired-end-bwa /
	                    paired-end-bfast)
	  -o out            output prefix

	Source FASTQ files should satisfy the following conditions: 1) Each file
	contains only reads corresponding to one genome (with the same genome id). 2)
	All files contain reads of the same type (single-end / paired-end). 3) Reads
	with more reads per tuple (e.g., paired-end) have '/1', etc. in suffix (for
	identification of nb of read).

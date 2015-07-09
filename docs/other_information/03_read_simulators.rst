.. _list-of-read-simulators:


Exhaustive list of read simulators
==================================

Here we provide an exhaustive list of existing read simulators in alphabetical order with additional information. Only several of them are supported by :ref:`MIShmash`.

.. contents::
	:depth: 3


Art
~~~


+--------------+-------------------------------------------------------------------------+
| Author:      | Weichun Huang                                                           |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://www.niehs.nih.gov/research/resources/software/biostatistics/art/ |
+--------------+-------------------------------------------------------------------------+
| Publication: | Huang, W. *et al.*                                                      |
|              | ART: a next-generation sequencing read simulator.                       |
|              | *Bioinformatics* **28**\(4), pp. 593--594, 2011.                        |
+--------------+-------------------------------------------------------------------------+



CuReSim
~~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Ségolène Caboche                                                        |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://www.pegase-biosciences.com/tools/curesim/                        |
+--------------+-------------------------------------------------------------------------+
| Publication: | Caboche, S. *et al.*                                                    |
|              | Comparison of mapping algorithms used in high-throughput sequencing:    |
|              | application to Ion Torrent data.                                        |
|              | *BMC Genomics* **15**\:264, 2014.                                       |
+--------------+-------------------------------------------------------------------------+


**Format of read names**

	.. code-block:: text

		@(.*)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)


1) contig name
2) original position
3) strand (0=forward;1=reverse)
4) random read (0=non-random;1=random)
5) number of insertions
6) number of deletions
7) number of substitution
8) read number (unique within a genome)


DNemulator
~~~~~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Martin C. Frith                                                         |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://cbrc3.cbrc.jp/~martin/dnemulator/                                |
+--------------+-------------------------------------------------------------------------+
| Publication: | Frith, M. C. *et al.*                                                   |
|              | A mostly traditional approach improves alignment of bisulfite-converted |
|              | DNA                                                                     |
|              | *Nucleic Acids Research* **40**\(13): e100, 2012.                       |
+--------------+-------------------------------------------------------------------------+




DwgSim
~~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Niels Homer                                                             |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://github.com/nh13/dwgsim                                           |
+--------------+-------------------------------------------------------------------------+

**Format of read names**

	.. code-block:: text

		(.*)_([0-9]+)_([0-9]+)_([01])_([01])_([01])_([01])_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_(([0-9abcdef])+)


1)  contig name (chromsome name)
2)  start end 1 (one-based)
3)  start end 2 (one-based)
4)  strand end 1 (0 - forward, 1 - reverse)
5)  strand end 2 (0 - forward, 1 - reverse)
6)  random read end 1 (0 - from the mutated reference, 1 - random)
7)  random read end 2 (0 - from the mutated reference, 1 - random)
8)  number of sequencing errors end 1 (color errors for colorspace)
9)  number of SNPs end 1
10) number of indels end 1
11) number of sequencing errors end 2 (color errors for colorspace)
12) number of SNPs end 2
13) number of indels end 2
14) read number (unique within a given contig/chromosome)



FastqSim
~~~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Anna Shcherbina                                                         |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://sourceforge.net/projects/fastqsim                                |
+--------------+-------------------------------------------------------------------------+
| Publication: | Shcherbina, A.                                                          |
|              | FASTQSim: platform-independent data characterization and in silico      |
|              | read generation for NGS datasets.                                       |
|              | *BMC Research Notes* **7**\:533, 2014.                                  |
+--------------+-------------------------------------------------------------------------+



FlowSim
~~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Susanne Balzer                                                          |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://biohaskell.org/Applications/FlowSim                              |
+--------------+-------------------------------------------------------------------------+
| Publication: | Balzer, S. *et al.*                                                     |
|              | Characteristics of 454 pyrosequencing data -- enabling realistic        |
|              | simulation with flowsim.                                                |
|              | *Bioinformatics* **26**\(18):i420--i425, 2010.                          |
+--------------+-------------------------------------------------------------------------+



GemSim
~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Kerensa E. McElroy                                                      |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://sourceforge.net/projects/gemsim                                  |
+--------------+-------------------------------------------------------------------------+
| Publication: | McElroy, K. E. *et al.*                                                 |
|              | GemSIM: general, error-model based simulator of next-generation         |
|              | sequencing data.                                                        |
|              | *BMC Genomics* **13**\:74, 2012.                                        |
+--------------+-------------------------------------------------------------------------+



Mason
~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Manuel Holtgrewe                                                        |
+--------------+-------------------------------------------------------------------------+
| URL:         | https://www.seqan.de/projects/mason/                                    |
+--------------+-------------------------------------------------------------------------+
| Publication: | Holtgrewe, M.                                                           |
|              | Mason -- a read simulator for second generation sequencing data.        |
|              | Technical Report TR-B-10-06,                                            |
|              | Institut für Mathematik und Informatik, Freie Universität Berlin, 2010. |
+--------------+-------------------------------------------------------------------------+



PbSim
~~~~~

+--------------+-------------------------------------------------------------------------+
| Authors:     | Michiaki Hamada, Yukiteru Ono                                           |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://code.google.com/p/pbsim                                          |
+--------------+-------------------------------------------------------------------------+
| Publication: | Ono, Y. *et al.*                                                        |
|              | PBSIM: PacBio reads simulator -- toward accurate genome assembly.       |
|              | *Bioinformatics* **29**\(1):119--121, 2013.                             |
+--------------+-------------------------------------------------------------------------+



Pirs
~~~~

+--------------+-------------------------------------------------------------------------+
| Authors:     | Jianying Yuan, Eric Biggers                                             |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://code.google.com/p/pirs                                           |
|              | http://github.com/galaxy001/pirs                                        |
+--------------+-------------------------------------------------------------------------+
| Publication: | Hu, X., *et al.*                                                        |
|              | pIRS: Profile-based Illumina pair-end reads simulator.                  |
|              | *Bioinformatics* **28**\(11):1533--1535, 2012.                          |
+--------------+-------------------------------------------------------------------------+



Sherman
~~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Felix Krueger                                                           |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://www.bioinformatics.babraham.ac.uk/projects/sherman               |
+--------------+-------------------------------------------------------------------------+



SimNgs
~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Tim Massingham                                                          |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://www.ebi.ac.uk/goldman-srv/simNGS/                                |
+--------------+-------------------------------------------------------------------------+



SimSeq
~~~~~~

+--------------+-------------------------------------------------------------------------+
| Authors:     | John St. John                                                           |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://github.com/jstjohn/SimSeq                                        |
+--------------+-------------------------------------------------------------------------+



SInC
~~~~

+--------------+-------------------------------------------------------------------------+
| URL:         | http://sincsimulator.sourceforge.net                                    |
+--------------+-------------------------------------------------------------------------+
| Publication: | Pattnaik, S. *et al.*                                                   |
|              | SInC: an accurate and fast error-model based simulator for SNPs, Indels |
|              | and CNVs coupled with a read generator for short-read sequence data.    |
|              | *BMC Bioinformatics* **15**\:40, 2014.                                  |
+--------------+-------------------------------------------------------------------------+


Wessim
~~~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Sangwoo Kim                                                             |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://sak042.github.io/Wessim/                                         |
+--------------+-------------------------------------------------------------------------+


 
WgSim
~~~~~

+--------------+-------------------------------------------------------------------------+
| Author:      | Heng Li                                                                 |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://github.com/lh3/wgsim                                             |
+--------------+-------------------------------------------------------------------------+


**Format of read names**

	.. code-block:: text

		(.*)_([0-9]+)_([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9a-f]+)/([12])



1) contig name (chromosome name)
2) start end 1 (one-based)
3) end end 2 (one-based)
4) number of errors end 1
5) number of substitutions end 1
6) number of indels end 1
7) number of errors end 2
8) number of substitutions end 2
9) number of indels end 2
10) id
11) pair



XS
~~

+--------------+-------------------------------------------------------------------------+
| Authors:     | Diogo Pratas, Armando J. Pinho, João M. O. S. Rodrigues                 |
+--------------+-------------------------------------------------------------------------+
| URL:         | http://bioinformatics.ua.pt/software/xs                                 |
+--------------+-------------------------------------------------------------------------+
| Publication: | Pratas, D. *et al.*                                                     |
|              | XS: a FASTQ read simulator.                                             |
|              | *BMC Research Notes* **7**\:40, 2014.                                   |
+--------------+-------------------------------------------------------------------------+

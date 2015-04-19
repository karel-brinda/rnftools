.. _introduction:

Introduction
============

RNFtools is an associate software package for `RNF`_, a generic format for assigning read names of simulated
Next-Generation Sequencing reads. The format aims to remove dependency of evaluation tools of read mappers
on the used read simulators. RNFtools consist of two principal parts:

* **MIShmash** - a tool for simulating NGS reads in RNF format using existing simulators.
* **LAVEnder** - a tool for evaluation of NGS read mappers using simulated reads in RNF format.

Technically, the entire RNFtools package is based on `SnakeMake`_, a `Make`_-like `Python`_-based software
developed primarly for creating bioinformatics pipelines. Basic usage of RNFtools is easily
comprehensible from `examples`_. For advance pipelines (e.g., simulating NGS reads + mapping + NGS mapper
evaluation in a single pipeline), some basic knowledge of Python is recommended.

RNFtools have few unique features:

* All required software and data files are installed fully automatically (using the `SMBL`_ library).
* Created pipelines are fully reproducible and they can be distributed as single SnakeMake files.


How to start with RNFtools
--------------------------

First we recommend to start with `tutorial`_ which demonstrates how to install RNFtools and how to use it. If anything is not clear, please look into `FAQ`_ if the answer already exists. If not, use the RNFtools mailing list. 


Have you found a bug?
---------------------

Please, report it as a `GitHub ticket`_. We will try to correct it as soon as possible.


.. _RNF: http://github.com/karel-brinda/rnf-spec/
.. _SMBL: http://github.com/karel-brinda/smbl/
.. _examples: http://github.com/karel-brinda/rnftools/tree/master/examples
.. _GitHub ticket: http://github.com/karel-brinda/rnftools/issues
.. _SnakeMake: http://bitbucket.org/johanneskoester/snakemake
.. _Make: http://www.gnu.org/software/make
.. _Python: http://python.org

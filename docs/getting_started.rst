.. _introduction:

Introduction
============

RNFtools is an associate software package for `RNF`_, a generic format for assigning read names of simulated
Next-Generation Sequencing reads. The format aims to remove dependency of evaluation tools of read mappers
on the used read simulators.

RNFtools consist of three principal parts:

1. :ref:`MIShmash` (MIS => SIM => simulation) - a tool for simulating NGS reads in RNF format using existing simulators.
2. :ref:`LAVEnder` (LAVE => EVAL => evaluation) - a tool for evaluation of NGS read mappers using simulated reads in RNF format.
3. :ref:`RNF` - a Python library for handling the RNF format.

Technically, the entire RNFtools package is based on `SnakeMake`_, a `Make`_-like `Python`_-based software developed for creating bioinformatics pipelines. It has few unique features:

	* All required software and data files are installed **fully automatically** (using the `SMBL`_ library).
	* Created pipelines are fully reproducible and they can be distributed as **single SnakeMake files**.

There exists also a console variant of RNFtools with a basic functionality (see :ref:`command_line`).

How to start with RNFtools
--------------------------

First we recommend to start with :ref:`tutorial` which demonstrates how to install RNFtools and how to use it. All examples are also located in a `dedicated directory`_. When you get familiar with basics, the main source of information will be :ref:`Reference`.

If anything is not clear, please look into :ref:`faq` and if this question has not been answered yet, send an e-mail to the `Mailing list`_.


.. _RNF: http://github.com/karel-brinda/rnf-spec/
.. _SMBL: http://github.com/karel-brinda/smbl/
.. _dedicated directory: http://github.com/karel-brinda/rnftools/tree/master/examples/tutorial
.. _GitHub ticket: http://github.com/karel-brinda/rnftools/issues
.. _SnakeMake: http://bitbucket.org/johanneskoester/snakemake
.. _Mailing list: http://groups.google.com/group/rnftools
.. _Make: http://www.gnu.org/software/make
.. _Python: http://python.org

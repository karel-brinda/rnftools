.. _extending_rnf:


Extending RNFtools
==================

This chapter devoted to developers describes how to do modifications in RNFtools.

A new simulator
---------------

To add support for another read simulator, you have to make the following four changes in source codes of RNFtools and SMBL. 


Step 1: SMBL plugin for installation of your simulator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First you should write a plugin for SMBL (https://github.com/karel-brinda/smbl) for automatic installation
of the simulator (plugins are located in this directory: https://github.com/karel-brinda/smbl/blob/master/smbl/prog/plugins/). You should get inspired from existing plugins such as https://github.com/karel-brinda/smbl/blob/master/smbl/prog/plugins/wgsim.py.
When the plugin is written, it should be added to https://github.com/karel-brinda/smbl/blob/master/smbl/prog/plugins/__init__.py.


Step 2: Wrapper of the simulator in RNFtools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All wrappers of read simulators are located in this directory: https://github.com/karel-brinda/rnftools/blob/master/rnftools/mishmash/. To create a
new wrapper, copy some existing one and modify the code inside.

Method ``__init__`` saves parameters for the simulation. Some of them are mandatory (read length, fasta file of the genome, etc.) and these are passed to ``__init__`` of the mother abstract class for a simulator.

Method ``get_input`` returns list of input files for simulation and one of them is your program (it should already exist in SMBL).

Method, ``create_fq`` simulates reads (by calling a shell command using ``smbl.utils.shell``) and converts the obtained files to RNF-FASTQ.

The previous step should be done using a static method (typically named ``recode_*_reads``) made for this purpose. All code for conversion should be located in this method such that it can be used externally without any instance of the class.

When all these functions are created/adjusted, the class should be imported in ``__init__.py`` located in the same directory. Get inspired by existing code, with a high probability only small changes will be needed.

At this point, it should be tested if the simulator works well in 


Step 3: Support in the ``rnftools`` program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The last step is pluging your converting static function into ``rnftools`` console program, which is contained in this file: https://github.com/karel-brinda/rnftools/blob/master/rnftools/scripts.py. You will  have to add a new subcommand (which will call the static function) and create a parser for it. Again, follow existing code.


Step 4: Tests
~~~~~~~~~~~~~

Add corresponding tests (see the test directory).



Remarks
~~~~~~~
Implicit support of RNF in a tool makes situation much easier. Then only the SMBL plugin must be created and a very simple RNFtools wrapper added. No conversion is needed.
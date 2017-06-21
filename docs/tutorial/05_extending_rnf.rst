.. _extending_rnf:

Extending RNFtools
==================

This chapter describes how to develop an extension for a new read simulator.


Step 1: Wrapper of the simulator in RNFtools
````````````````````````````````````````````

All wrappers of read simulators are located in this directory:
https://github.com/karel-brinda/rnftools/blob/master/rnftools/mishmash/. To
create a new wrapper, copy some existing one and modify the code inside.

Method ``__init__`` saves the parameters for a simulation. Some of them are
mandatory (read length, FASTA file of the genome, etc.) and these are passed to
``__init__`` of the mother abstract class for a simulator.

Method ``get_input`` returns list of input files for simulation and one of them
is your program.

Method, ``create_fq`` simulates reads (by calling a shell command using
``rnftools.utils.shell``) and converts the obtained files to RNF-FASTQ.

The previous step should be done using a dedicated static method (typically
named ``recode_*_reads``). All code for conversion should be located in this
method so that it can be used externally without creating an instance of the
class.

When all these functions are created/adjusted, the class should be imported in
``__init__.py``. Get inspired by existing code, with a high probability only
little changes will be needed.

.. At this point, it should be tested if the simulator works well in


Step 2: Support in the ``rnftools`` program
```````````````````````````````````````````

The last step is plugging your converting static function into ``rnftools``
console program, which is contained in this file:
https://github.com/karel-brinda/rnftools/blob/master/rnftools/scripts.py. You
will  have to add a new subcommand (which will call the static function) and
create a parser for it. Again, follow existing code.


Step 3: Tests
`````````````

Add corresponding tests (see the test directory).


Information for developers
==========================

Do you develop bioinformatics software? Here you will find how RNF and RNFtools can be useful for you .

.. contents::
   :depth: 3


Do you develop a **read mapper**?
---------------------------------

RNFtools can help you to debug your mapper. You can: 

* Find reads which were not aligned correctly (when tuning your algorithm).
* Test how successful your mapper is in dependence on parameters of simulation (error rate, etc.).
* Test if contamination reads are well detected and staying really unaligned.
* Test which impact have pre- and post-processing tools (such as read clustering tools or re-alignment tools) in combination with you mapper.

First you can start with some simple simulator (e.g., WgSim) for basic tests and later easily switch to more realistic simulations with ART of Mason.


Do you develop a **read simulator**?
------------------------------------

Even though MIShmash currently supports several simulators, implicit support of RNF in simulators is preferable. Usually, simulators do not save as much information as RNF does, hence MIShmash must sometimes estimate some of these unknown values which brings noise into data. Direct support of RNF support in simulators would imply higher precision in the forthcoming analysis as well as better usability of your software.

Adding support for RNF into your simulator is a simple step since the format is easy to adopt. For existing software, we recommend to add an extra parameter switching internal naming procedure to RNF.

RNFtools can also help you to debug your simulator. Switching between your simulator and another one, you can check if obtained results are similar as they should be. If not, you might have bugs in your code (such as coordinates of simulation are incorrectly shifted by few positions, etc.).


Do you develop an **evaluation tool for read mappers**?
-------------------------------------------------------

RNF enables you writing a universal evaluation tool, compatible with all RNF-compatible read simulators and all simulators supported by MIShmash.

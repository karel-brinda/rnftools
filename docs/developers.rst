Information for developers
==========================

Do you develop bioinformatics software? Here you will find how RNF and RNFtools
can be useful.


\... of mappers of NGS reads
----------------------------

RNFtools can help you to debug your tool. You can simulate reads using MIShmash,
map them using your mapper, and then evaluate it using LAVEnder.

Then you can modify parameters for the simulator and also for your mapper and
observe how these changes affect sensitivity and precision of mapping.

Then you can also easily swithc between the supported read simulators. It can be
contributive since different simulators can simulate different artefacts in data.

In similar way, RNFtools can be used for testing of pre-processing and
post-processing NGS tools (such as read clustering tools or re-alignment tools).


\... of simulators of NGS reads
-------------------------------

Even though MIShmash currently supports several simulators, it is much better when
a tool supports the RNF naming convention implicitely. Usually, simulators do not
save as much information as it is saved by RNF, hence MIShmash must sometimes estimate
some of these values values. Therefore, direct support in a read simulator implies higher
precision in the forthcoming analysis. It also increases usability of your software.

Adding support for RNF into your simulator is a simple step since the format is easy
to adopt. For existing software, it is suggested to add an extra parameter which
switches from the default naming convention to RNF.

RNFtools can also help you to debug your simulator. You can easily switch from your
simulator to another and compare results. When the parameters are similar, also the 
obtained results should be similar. You can also find many problems using
mappers -- by comparing if the mappers really do map the reads to the expected
positions. E.g., many reads are shifted by two positions, there is probably a bug
inside.


\... of evaluators of read mappers
----------------------------------

RNF enables you writing a universal evaluation tool, compatible with any RNF-compatible read
simulator and any simulator supported by MIShmash.

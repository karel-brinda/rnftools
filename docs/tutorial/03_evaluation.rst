Mapper evalution
----------------

In this chapter, we show how to evaluate read mappers using RNFtools.
For this task, we will use a component called LAVEnder.

Basic example
"""""""""""""

The basic approach of mapper evaluation consists of the following steps:

1. Simulation of reads (see the previous chapter of this tutorial).
2. Mapping reads to a reference genome.
3. Creating the report.

The first step was described in the previous chapter. Once your simulated reads
are mapped (step 2), you can create the following Snakefile:

.. literalinclude:: ../../examples/01_tutorial/03_evaluation/03_evaluation/Snakefile
	:language: python
	:linenos:

When you run `snakemake`, RNFtools detect all BAM files in the specified
directories, and starts evaluation and creates an interactive HTML report
containing one panel for each directory.  The ``name`` argument defines the
name of the report (the final HTML file will have name ``{name}.html``).
Parameter ``keep_intermediate_files`` sets if the intermediate ES (evaluated
individual segments) and ET (evaluated read pairs) files created during
evaluation should be kept. The argument ``allowed_delta`` is used for setting
maximum allowed distance between reported position and original position for
considering the segment still correctly mapped.


Auxiliary files
"""""""""""""""

For every BAM file, the following files are created.

* **HTML** -- detailed report for the BAM file
* **ES** (mapping information: segments) -- file with information about mapping categories of each segment
* **ET** (mapping information: read tuples) -- file with information about category of entire read tuples
* **ROC** -- source file for plotting graphs, it contains information about how many reads are in which category in dependence on threshold on mapping qualities
* **GP** -- GnuPlot file used for plotting detailed graphs for this BAM  (**SVG**, **PDF**)


Adjusting plotted graphs
""""""""""""""""""""""""

For details about adjusting graphs, please see :ref:`lavender`.

First, you can change default values for the basic graphs:

* range of x and y axes (``default_x_run``, ``default_y_run``),
* sizes of created PDF and SVG files (``default_pdf_size_cm``, ``default_svg_size_px``)
* label of x-axis (``default_x_label``)
* values on x-axis (``default_x_axis``).

Then you can add your own graphs using ``rnftools.lavender.Report.add_graph`` method.



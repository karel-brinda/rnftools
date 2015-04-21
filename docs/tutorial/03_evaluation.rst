Mapper evalution
----------------

In this chapter we show how to evaluate read mappers using read names in RNF format. For this task, component called LAVEnder is used.

Basic example
"""""""""""""

The basic approach of mapper evaluation consists of the following steps:

1. Simulation of reads (see the previous chapter of this tutorial).
2. Mapping reads to a reference genome.
3. Creating the report.

First you need to simulate RNF reads and map them to a reference genome. If you don't have any own BAM file, please use the following toy Snakefile which performs first two steps:

.. literalinclude:: ../../examples/tutorial/03_evaluation/01_simple_evaluation/bams/Snakefile
	:language: python
	:linenos:

In a directory for this experiment, create a directory ``bams`` and place there the previous code and run ``snakemake`` there. If you have your own BAM files, create the ``bams`` directory as well and place them there.

Now let us create a directory ``report`` with the following ``Snakefile``:

.. literalinclude:: ../../examples/tutorial/03_evaluation/01_simple_evaluation/report/Snakefile
	:language: python
	:linenos:

Evaluation of all BAM files in a dir is requested by creating an instance of the ``rnftools.lavender.Report`` class (line 3). Parameter ``bam_dirs`` accepts a list of directories with BAM files. Every entry of the list corresponds to a single panel in the final HTML report. The ``name`` argument defines name of the report (the final HTML file will have name ``{name}.html``). Parameter ``keep_intermediate_files`` indicates if intermediate MIS and MIR files created during evaluation should be kept. Argument ``allowed_delta`` is used for setting maximum allowed distance between reported position and original position for considering the segment still correctly mapped. 

When you execute ``snakemake``, the report is created.


Auxiliary files
"""""""""""""""

For every BAM file, the following files are created.

* **HTML** -- detailed report for the BAM file
* **MIS** (mapping information: segments) -- file with information about mapping categories of each segment
* **MIR** (mapping information: read tuples) -- file with information about category of entire read tuples
* **ROC** -- source file for plotting graphs, it contains information about how many reads are in which category in dependence on threshold on mapping qualities
* **GP** -- GnuPlot file used for plotting detailed graphs for this BAM  (**SVG**, **PDF**)


Adjusting plotted graphs
""""""""""""""""""""""""

For details about adjusting graphs, please see Reference/LAVEnder.

First, you can change default values for the basic graphs:

* range of x and y axes (``default_x_run``, ``default_y_run``), 
* sizes of created PDF and SVG files (``default_pdf_size_cm``, ``default_svg_size_px``)
* label of x-axis (``default_x_label``)
* values on x-axis (``default_x_axis``).

Then you can add your own graphs using ``rnftools.lavender.Report.add_graph`` method.



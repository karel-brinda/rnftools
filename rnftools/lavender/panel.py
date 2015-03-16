import rnftools.lavender

import smbl
import snakemake
import os
import glob

__all__=["Panel"]

#############
#############
### PANEL ###
#############
#############
class Panel:
	"""Class for a single panel in a HTML report."""

	def __init__(self,
		report,
		bam_dir,
		panel_dir,
		name
		):
		"""

		:param report: The owner (report).
		:type  report: str
		:param bam_dir: Directory to the BAM files.
		:type  bam_dir: str
		:param panel_dir: Directory with panel's auxiliary files.
		:type  panel_dir: str
		:param name: Name of the panel.
		:type  name: str
		:raises: ValueError

		"""
		self.report=report
		rnftools.lavender.add_panel(self)
		self.name=name
		self.panel_dir=panel_dir

		self._svg_fn=os.path.join(self.panel_dir,"svg","_combined.svg")
		self._gp_fn=os.path.join(self.panel_dir,"gp","_combined.gp")
		self._pdf_fn=os.path.join(self.panel_dir,"pdf","_combined.pdf")

		bams_fns=glob.glob(os.path.join(bam_dir,"*.bam"))
		self.bams=[
				rnftools.lavender.Bam(
					bam_fn=bam_fn,
					panel=self,
					name=os.path.basename(bam_fn).replace(".bam","")
				) 
				for bam_fn in sorted(bams_fns)
			]

		if len(self.bams)==0:
			raise ValueError("Panel '{}' does not contain any BAM file.".format(self.name))

		for x in ["gp","html","roc","svg"]:
			snakemake.shell('mkdir -p "{}"'.format(os.path.join(self.panel_dir,x)))

	def get_report(self):
		""" Get the report. """

		return self.report

	def get_panel_dir(self):
		""" Get the directory with panel's auxiliary files. """

		return self.panel_dir

	def get_bams(self):
		""" Get BAMs for this panel. """

		return self.bams

	def get_required_files(self):
		""" Get all required files. """

		return [bam.get_required_files() for bam in self.bams] + [self._svg_fn]

	def get_html_column(self):
		""" Get a HTML column for this panel. """

		panel_id="panel_{}".format(self.name)
		return [
				(" <br>"+os.linesep).join(
					[
						"""
							<strong>{bam_name}:</strong>
							<a onclick="document.getElementById('{panel_id}').src='{bam_svg}';return false;" href="#">graph</a>,
							<a href="{bam_html}">report</a>:
						""".format(
							bam_name=bam.get_name(),
							bam_html=bam.html_fn(),
							bam_svg=bam.svg_fn(),
							panel_id=panel_id,
						)
						for bam in self.bams
					]

				),

				"""<img src="{svg}" id="{panel_id}">""".format(
						svg=self.bams[0]._svg_fn,
						panel_id=panel_id
					),

				"""<img src="{svg}">""".format(
						svg=self._svg_fn
					)
			]

	######################################
	######################################

	def gp_fn(self):
		""" Get the GnuPlot file name for the overall graph. """

		return self._gp_fn

	def svg_fn(self):
		""" Get the SVG file name for the overall graph. """

		return self._svg_fn

	def pdf_fn(self):
		""" Get the PDF file name for the overall graph. """

		return self._pdf_fn

	######################################
	######################################

	def create_gp(self):
		""" Create GnuPlot file. """

		# todo: parameter
		def gp_style(i):
			colors=["red","green","blue","goldenrod","black"]
			color=colors[i % len(colors)]
			return 'set style line {i} lt 1 pt {i} lc rgb "{color}";'.format(color=color,i=i+1)

		with open(self._gp_fn,"w+") as f:
			plots =[
				""""{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2+$3+$4)*100/($2+$3+$4+$7+$8)) \
					with lp ls {i} ps 0.8 title "{basename}" noenhanced,\\""".format(
							roc_fn=self.bams[i].roc_fn(),
							basename=os.path.basename(self.bams[i].roc_fn()),
							i=i+1,
						)
					for i in range(len(self.bams))
			]

			f.write("""
				set key spacing 0.8 opaque width -3

				set x2lab "false positive rate\\n(#wrong mappings / #mapped)"
				set log x
				set log x2

				set format x "10^{{%L}}"
				set format x2 "10^{{%L}}"
				set xran  [{xran}]
				set x2ran [{xran}]
				set x2tics

				{styles}

				set ylab "sensitivity on reads to map (%)"

				set format y "%g %%"
				set yran [60:100]

				set pointsize 1.5

				set grid ytics lc rgb "#777777" lw 1 lt 0 front
				set grid xtics lc rgb "#777777" lw 1 lt 0 front

				set datafile separator "\t"
				set palette negative

				set termin svg size {svg_size} enhanced
				set out "{svg_fn}"
				{plots}


				set termin pdf enhanced size {pdf_size} enhanced font 'Arial,12'
				set out "{pdf_fn}"
				{plots}

				""".format(
					svg_fn=self._svg_fn,
					pdf_fn=self._pdf_fn,
					xran="{:.10f}:{:.10f}".format(self.report.plot_x_run[0],self.report.plot_x_run[1]),
					yran="{:.10f}:{:.10f}".format(self.report.plot_y_run[0],self.report.plot_y_run[1]),
					svg_size="{},{}".format(self.report.plot_svg_size[0],self.report.plot_svg_size[1]),
					pdf_size="{:.10f}cm,{:.10f}cm".format(self.report.plot_pdf_size_cm[0],self.report.plot_pdf_size_cm[1]),
					plots="plot "+os.linesep.join(plots),
					styles=os.linesep.join([gp_style(i) for i in range(40)]),
				)
			)

	def create_graphics(self):
		"""Create images related to this panel."""

		snakemake.shell('"{}" "{}"'.format(smbl.prog.GNUPLOT5,self._gp_fn))

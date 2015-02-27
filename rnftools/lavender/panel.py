import lavender
import smbl
import snakemake
import os
import glob

##############
##############
### PANNEL ###
##############
##############
class Pannel:
	"""Class for a single pannel in a HTML report."""

	def __init__(self,report,bam_dir,pannel_dir,name):
		"""

		:param report: The owner (report).
		:param bam_dir: Directory to the BAM files.
		:param pannel_dir: Directory with pannel's auxiliary files.
		:param name: Name of the pannel.

		"""
		self.report=report
		lavender._PANNELS_.append(self)
		self.name       = name
		self.pannel_dir = pannel_dir

		self._svg_fn    = os.path.join(self.pannel_dir,"svg","_combined.svg")
		self._gp_fn     = os.path.join(self.pannel_dir,"gp","_combined.gp")

		bams_fns        = glob.glob(os.path.join(bam_dir,"*.bam"))
		self.bams       = [
				lavender.Bam(
					bam_fn=bam_fn,
					pannel=self,
					name=os.path.basename(bam_fn).replace(".bam","")
				) 
				for bam_fn in sorted(bams_fns)
			]

		for x in ["gp","html","roc","svg"]:
			snakemake.shell('mkdir -p "{}"'.format(os.path.join(self.pannel_dir,x)))

	def get_report(self):
		""" Get the report. """

		return self.report

	def get_pannel_dir(self):
		""" Get the directory with pannel's auxiliary files. """

		return self.pannel_dir

	def get_bams(self):
		""" Get BAMs for this pannel. """

		return self.bams

	def get_required_files(self):
		""" Get all required files. """

		return [bam.get_required_files() for bam in self.bams] + [self._svg_fn]

	def get_html_column(self):
		""" Get a HTML column for this pannel. """

		pannel_id="pannel_{}".format(self.name)
		return [
				(" <br>"+os.linesep).join(
					[
						"""
							<a style="font-weight:bold" href="{bam_html}">{bam_name}</a>: 
							<a onclick="document.getElementById('{pannel_id}').src='{bam_svg}';return false;" href="#">graph</a>
						""".format(
							bam_name=bam.get_name(),
							bam_html=bam.html_fn(),
							bam_svg=bam.svg_fn(),
							pannel_id=pannel_id
						)
						for bam in self.bams
					]

				),

				"""<img src="{svg}" id="{pannel_id}">""".format(
						svg=self.bams[0]._svg_fn,
						pannel_id=pannel_id
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

	######################################
	######################################

	def create_gp(self):
		""" Create GnuPlot file. """

		with open(self._gp_fn,"w+") as f:
			plots =[
				""""{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2)*100/$10) \
					with line dt 1 lc rgb "#0000ff" title "{basename}" noenhanced,\\""".format(
							roc_fn=bam.roc_fn(),
							basename=os.path.basename(bam.roc_fn())
						)
					for bam in self.bams
			]

			f.write("""
				set ylab "part of all reads (%)"

				set x2lab "1 - precision\\n(#wrong mappings / #mapped)"
				set log x
				set log x2
				set xran [0.00001:1]
				set x2ran [0.00001:1]
				set x2tics

				set yran [60:100]
				set format y "%g %%"

				set pointsize 1.5

				set termin svg size 640,640
				set out "{svg}"

				set grid ytics lc rgb "#777777" lw 1 lt 0 front
				set grid xtics lc rgb "#777777" lw 1 lt 0 front

				set datafile separator "\t"

				plot \\
            	{plots}
				""".format(
					plots=os.linesep.join(plots),
					svg=self._svg_fn
				)
			)

	def create_svg(self):
		""" Create SVG file. """

		snakemake.shell('{} "{}"'.format(smbl.prog.GNUPLOT5,self._gp_fn))

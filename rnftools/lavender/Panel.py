import rnftools.lavender

import smbl
import snakemake
import os
import glob

from . import _default_gp_style_func

#############
#############
### PANEL ###
#############
#############
class Panel:
	"""Class for a single panel in a HTML report.

	Args:
		report (rnftools.lavender.Report): The owner report.
		bam_dir (str): Directory to the BAM files for this panel.
		panel_dir (str): Directory with auxiliary files for this panel.
		name (str): Name of the panel.
		keep_intermediate_files (bool): Keep files created in intermediate steps during evaluation.
		compress_intermediate_files (bool): Compress files created in intermediate steps during evaluation.
		default_x_axis (str): Values on x-axis, e.g., "({m}+{w})/({M}+{m}+{w})".
		default_x_label (str): Label on x-axis.
		gp_style_func (function(i, nb)): Function assigning GnuPlot styles for overall graphs. Arguments: i: 0-based id of curve, nb: number of curves.
	
	Raises:
		ValueError

	"""

	def __init__(self,
			report,
			bam_dir,
			panel_dir,
			name,
			keep_intermediate_files,
			compress_intermediate_files,
			default_x_axis,
			default_x_label,
			gp_style_func=_default_gp_style_func,
		):

		self.report=report
		rnftools.lavender.add_panel(self)
		self.name=name
		self.panel_dir=panel_dir
		self.default_x_axis=default_x_axis
		self.default_x_label=default_x_label

		self.gp_plots = []

		self._gp_style_func=gp_style_func
		#simple test
		nb_styles=10
		for i in range(nb_styles):
			res=self._gp_style_func(i,nb_styles)
			assert isinstance(res,str)
			assert len(res)>0

		self._gp_fn=os.path.join(self.panel_dir,"gp","_combined.gp")
		self._svg_fns=[] #os.path.join(self.panel_dir,"svg","_combined.svg")
		self._pdf_fns=[] #os.path.join(self.panel_dir,"pdf","_combined.pdf")

		bams_fns=glob.glob(os.path.join(bam_dir,"*.bam"))
		self.bams=[
				rnftools.lavender.Bam(
					bam_fn=bam_fn,
					panel=self,
					name=os.path.basename(bam_fn).replace(".bam",""),
					keep_intermediate_files=keep_intermediate_files,
					compress_intermediate_files=compress_intermediate_files,
					default_x_axis=default_x_axis,
					default_x_label=default_x_label,
				) 
				for bam_fn in sorted(bams_fns)
			]

		if len(self.bams)==0:
			raise ValueError("Panel '{}' does not contain any BAM file.".format(self.name))

		for x in ["gp","html","roc","svg","pdf"]:
			smbl.utils.shell('mkdir -p "{}"'.format(os.path.join(self.panel_dir,x)))

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

		return [bam.get_required_files() for bam in self.bams] + [self._svg_fns]

	def get_html_column(self):
		""" Get a HTML column for this panel. """

		panel_id="panel_{}".format(self.name)
		return [
				(" <br />"+os.linesep).join(
					[
						"""
							<strong>{bam_name}:</strong>
							<a onclick="document.getElementById('{panel_id}').src='{bam_svg}';document.getElementById('{panel_id}_').href='{bam_html}';return false;" href="#">display graph</a>,
							<a href="{bam_html}">detailed report</a>
						""".format(
							bam_name=bam.get_name(),
							bam_html=bam.html_fn(),
							bam_svg=bam.svg_fn(),
							panel_id=panel_id,
						)
						for bam in self.bams
					]

				),

				"""
					<div class="formats">
						<a href="{html}" id="{panel_id}_">
							<img src="{svg}" id="{panel_id}" />
						</a>
					</div>
				""".format(
						html=self.bams[0]._html_fn,
						svg=self.bams[0]._svg_fn,
						pdf=self.bams[0]._pdf_fn,
						panel_id=panel_id,
					),
			] + [
				"""
					<div class="formats">
						<img src="{svg}" />
						<br />
						<a href="{pdf}">PDF version</a>
						|
						<a href="{svg}">SVG version</a>
						|
						<a href="{gp}" type="text/plain">GP file</a>
					</div>

				""".format(
						svg=svg,
						pdf=pdf,
						gp=self._gp_fn,
					)

				for (svg,pdf) in zip(self._svg_fns,self._pdf_fns)
			]
			

	######################################
	######################################

	def gp_fn(self):
		""" Get the GnuPlot file name for the overall graphs. """

		return self._gp_fn

	def svg_fns(self):
		""" Get the SVG file names for the overall graphs. """

		return self._svg_fns

	def pdf_fns(self):
		""" Get the PDF file names for the overall graphs. """
		return self._pdf_fns

	######################################
	######################################

	def add_graph(self,
				y,
				y_label,
				x_label,
				title,
				x_run,
				y_run,
				pdf_size_cm,
				svg_size_px,
				key_position,
			):

		x_gp=rnftools.lavender._format_xxx(self.default_x_axis)
		y_gp=rnftools.lavender._format_xxx("({})*100".format(y))	

		i=len(self.gp_plots)
		svg_file = os.path.join(self.panel_dir,"svg","_combined_{}.svg".format(i))
		pdf_file = os.path.join(self.panel_dir,"pdf","_combined_{}.pdf".format(i))

		self._svg_fns.append(svg_file)
		self._pdf_fns.append(pdf_file)

		params = [
			"",
			'set title "{{/:Bold {}}}"'.format(title),
			'set key {}'.format(key_position),
			'set x2lab "{}"'.format(x_label),
			'set ylab "{}"'.format(y_label),
			'set xran [{xran}]'.format(
						xran="{:.10f}:{:.10f}".format(x_run[0],x_run[1])
					),
			'set x2ran [{xran}]'.format(
						xran="{:.10f}:{:.10f}".format(x_run[0],x_run[1])
					),
			'set yran [{yran}]'.format(					
						yran="{:.10f}:{:.10f}".format(y_run[0],y_run[1]),
					),
			'set y2ran [{yran}]'.format(					
						yran="{:.10f}:{:.10f}".format(y_run[0],y_run[1]),
					),
			"",
		]

		plot =	[
					""""{roc_fn}" using ({x}):({y}) \
						with lp ls {i} ps 0.8 title "  {basename}" noenhanced,\\""".format(
								x=x_gp,
								y=y_gp,
								roc_fn=self.bams[i].roc_fn(),
								basename=os.path.basename(self.bams[i].roc_fn())[:-4],
								i=i+1,
							)
						for i in range(len(self.bams))
			]

		self.gp_plots.append( os.linesep.join(
				[
					"set termin pdf enhanced size {pdf_size} enhanced font 'Arial,12'".format(
							pdf_size="{:.10f}cm,{:.10f}cm".format(pdf_size_cm[0],pdf_size_cm[1])
						),
					'set out "{}"'.format(pdf_file),
					'set key spacing 0.8 opaque',
				] + params + [

					"plot \\"
				] + plot + ["",""]
			))

		self.gp_plots.append( os.linesep.join(
				[
					"set termin svg size {svg_size} enhanced".format(
							svg_size="{},{}".format(svg_size_px[0],svg_size_px[1])
						),
					'set out "{}"'.format(svg_file),
					'set key spacing 0.8 opaque',
				] + params + [


					"plot \\"
				] + plot + ["",""]
			))


	def create_gp(self):
		""" Create GnuPlot file. """

		nb_bams=len(self.bams)

		with open(self._gp_fn,"w+") as f:

			f.write("""
				set log x
				set log x2


				#set format x "10^{{%L}}"
				set format x2 "10^{{%L}}"
				set x2tics
				unset xtics

				{styles}

				set format y "%g %%"
				set ytics

				set pointsize 1.5

				set grid ytics lc rgb "#777777" lw 1 lt 0 front
				set grid x2tics lc rgb "#777777" lw 1 lt 0 front

				set datafile separator "\t"
				set palette negative

				{all_plots}

				""".format(
					all_plots=os.linesep.join(self.gp_plots),
					styles=os.linesep.join([self._gp_style_func(i,nb_bams) for i in range(nb_bams)]),
					x_lab=self.default_x_label,
				)
			)

	def create_graphics(self):
		"""Create images related to this panel."""

		if len(self._svg_fns)+len(self._pdf_fns)>0:
			smbl.utils.shell('"{}" "{}"'.format(smbl.prog.GNUPLOT5,self._gp_fn))

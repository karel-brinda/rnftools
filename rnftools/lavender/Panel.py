import rnftools.lavender
import rnftools.utils

import snakemake
import os
import glob
import tarfile
import io
import re
import textwrap

from . import _default_gp_style_func
from svg42pdf import svg42pdf


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
		name (str): Name of the panel (used for CSS, etc.).
		title (str): Title of the panel (to be displayed).
		render_pdf_method (str): Method for svg42pdf to render PDF (None / 'any' / 'cairo' / 'reportlab' / 'inkscape' / 'imagemagick' / 'wkhtmltopdf').
		keep_intermediate_files (bool): Keep files created in intermediate steps during evaluation.
		compress_intermediate_files (bool): Compress files created in intermediate steps during evaluation.
		default_x_axis (str): Values on x-axis, e.g., "({m}+{w})/({M}+{m}+{w})".
		default_x_label (str): Label on x-axis.
		gp_style_func (function(i, nb)): Function assigning GnuPlot styles for overall graphs. Arguments: i: 0-based id of curve, nb: number of curves.

	Raises:
		ValueError

	"""

	def __init__(
			self,
			report,
			bam_dir,
			panel_dir,
			name,
			title,
			render_pdf_method,
			keep_intermediate_files,
			compress_intermediate_files,
			default_x_axis,
			default_x_label,
			gp_style_func=_default_gp_style_func,
	):

		self.report = report
		rnftools.lavender.add_panel(self)
		self.name = name
		self.title = title
		self.panel_dir = panel_dir
		self.default_x_axis = default_x_axis
		self.default_x_label = default_x_label

		self.render_pdf_method = render_pdf_method

		self.gp_plots = []

		self._gp_style_func = gp_style_func
		# simple test
		nb_styles = 10
		for i in range(nb_styles):
			res = self._gp_style_func(i, nb_styles)
			assert isinstance(res, str)
			assert len(res) > 0

		self._gp_fn = os.path.join(self.panel_dir, "gp", "_combined.gp")
		self._tar_fn = os.path.join(self.panel_dir, "tar",
			"{title}.{panel}.tar".format(title=self.report.title, panel=self.title))
		self._svg_fns = []
		self._pdf_fns = []

		bams_fns = glob.glob(os.path.join(bam_dir, "*.bam"))
		self.bams = [
			rnftools.lavender.Bam(
				bam_fn=bam_fn,
				panel=self,
				name=os.path.basename(bam_fn).replace(".bam", ""),
				render_pdf_method=self.render_pdf_method,
				keep_intermediate_files=keep_intermediate_files,
				compress_intermediate_files=compress_intermediate_files,
				default_x_axis=default_x_axis,
				default_x_label=default_x_label,
			)
			for bam_fn in sorted(bams_fns)
		]

		if len(self.bams) == 0:
			raise ValueError("Panel '{}' does not contain any BAM file.".format(self.name))

		for x in ["gp", "html", "roc", "graphics", "tar"]:
			rnftools.utils.shell('mkdir -p "{}"'.format(os.path.join(self.panel_dir, x)))

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

		panel_id = "panel_{}".format(self.name)
		return [
				   "<h2>{}</h2>".format(self.title) +
				   '<a href="{}">Download data</a>'.format(self.tar_fn())
			   ] + [
				   # list of links
				   (" <br />" + os.linesep).join(
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

				   ) + '<br /> '.format(self.tar_fn()),

				   # main graph
				   """
				   <div class="formats">
				   <a href="{html}" id="{panel_id}_">
				   <img src="{svg}" id="{panel_id}" />
				   </a>
				   </div>
				   """.format(
					   html=self.bams[0]._html_fn,
					   svg=self.bams[0]._svg_fn,
					   panel_id=panel_id,
				   ),
			   ] + [
				   # overall graphs
				   """
				   <div class="formats">
				   <img src="{svg}" />
				   <br />
				   <a href="{svg}">SVG version</a>
				   |
				   <a href="{gp}" type="text/plain">GP file</a>
				   </div>
				   
				   """.format(
					   svg=svg,
					   gp=self._gp_fn,
				   )

				   for svg in self._svg_fns
			   ]

	######################################
	######################################

	def tar_fn(self):
		""" Get the TAR file name for the archive with data. """

		return self._tar_fn

	def gp_fn(self):
		""" Get the GnuPlot file name for the overall graphs. """

		return self._gp_fn

	def svg_fns(self):
		""" Get the PDF file names for the overall graphs (empty list if they are not rendered). """

		return self._svg_fns

	def pdf_fns(self):
		if self.render_pdf_method is None:
			return []
		else:
			return [re.sub(r'\.svg$', r'.pdf', svg_fn) for svg_fn in self._svg_fns]

	######################################
	######################################

	def add_graph(
			self,
			y,
			y_label,
			x_label,
			title,
			x_run,
			y_run,
			svg_size_px,
			key_position,
	):

		x_gp = rnftools.lavender._format_xxx(self.default_x_axis)
		y_gp = rnftools.lavender._format_xxx("({})*100".format(y))

		i = len(self.gp_plots)
		svg_file = os.path.join(self.panel_dir, "graphics", "_combined_{}.svg".format(i))

		self._svg_fns.append(svg_file)

		params = [
			"##################################################",
			'set title "{{/:Bold {}}}"'.format(title),
			'set key {}'.format(key_position),
			'set x2lab "{}"'.format(x_label),
			'set ylab "{}"'.format(y_label),
			'set xran [{xran}]'.format(
				xran="{:.10f}:{:.10f}".format(x_run[0], x_run[1])
			),
			'set x2ran [{xran}]'.format(
				xran="{:.10f}:{:.10f}".format(x_run[0], x_run[1])
			),
			'set yran [{yran}]'.format(
				yran="{:.10f}:{:.10f}".format(y_run[0], y_run[1]),
			),
			'set y2ran [{yran}]'.format(
				yran="{:.10f}:{:.10f}".format(y_run[0], y_run[1]),
			),
			"",
		]

		plot = [
			""""{roc_fn}" using ({x}):({y}) with lp ls {i} ps 0.8 title "  {basename}" noenhanced,\\""".format(
				x=x_gp,
				y=y_gp,
				roc_fn=self.bams[i].roc_fn(),
				basename=os.path.basename(self.bams[i].roc_fn())[:-4],
				i=i + 1,
			)
			for i in range(len(self.bams))
		]

		self.gp_plots.append(os.linesep.join(
			[
				"set termin svg size {svg_size} enhanced".format(
					svg_size="{},{}".format(svg_size_px[0], svg_size_px[1])
				),
				'set out "{}"'.format(svg_file),
				'set key spacing 0.8 opaque',
			] + params + [

				"plot \\"
			] + plot + ["", ""]
		))

	def create_gp(self):
		""" Create GnuPlot file. """

		nb_bams = len(self.bams)

		gp_parts = [
			textwrap.dedent("""\
				set log x
				set log x2


				#set format x "10^{{%L}}"
				set format x2 "10^{{%L}}"
				set x2tics
				unset xtics
				"""
			),

			os.linesep.join([self._gp_style_func(i, nb_bams) for i in range(nb_bams)]),

			textwrap.dedent("""\
					set format y "%g %%"
					set ytics

					set pointsize 1.5

					set grid ytics lc rgb "#777777" lw 1 lt 0 front
					set grid x2tics lc rgb "#777777" lw 1 lt 0 front

					set datafile separator "\\t"
					set palette negative
					"""
			),

			os.linesep.join(self.gp_plots)
		]

		gp_src = os.linesep.join(gp_parts)
		# .format(
		# 	x_lab=self.default_x_label,
		# )

		with open(self._gp_fn, "w+") as f:
			f.write(gp_src)

	def create_graphics(self):
		"""Create images related to this panel."""

		if len(self._svg_fns) > 0:
			rnftools.utils.shell('"{}" "{}"'.format("gnuplot", self._gp_fn))

			if self.render_pdf_method is not None:
				for svg_fn in self._svg_fns:
					pdf_fn = re.sub(r'\.svg$', r'.pdf', svg_fn)
					svg42pdf(svg_fn, pdf_fn, method=self.render_pdf_method)

	def create_tar(self):
		"""Create a tar file with all the files."""

		def add_file_to_tar(tar, orig_fn, new_fn, func=None):
			tf = tarfile.TarInfo(name=new_fn)
			with open(orig_fn) as f:
				tfs = f.read()

			if func is not None:
				tfs = func(tfs)
			tf.size = len(tfs)
			tfs = io.BytesIO(tfs.encode('utf8'))
			tar.addfile(tarinfo=tf, fileobj=tfs)

		def add_text_to_tar(tar, new_fn, text, func=None):
			tf = tarfile.TarInfo(name=new_fn)
			if func is not None:
				text = func(text)
			tf.size = len(text)
			tfs = io.BytesIO(text.encode('utf8'))
			tar.addfile(tarinfo=tf, fileobj=tfs)

		def strip_lines(text):
			text = text.replace("\t", " ")
			while text.find("  ") != -1:
				text = text.replace("  ", " ")
			lines = [x.strip() for x in text.strip().split("\n")]
			return "\n".join(lines) + "\n"

		tar = tarfile.TarFile(self._tar_fn, "w")
		for i in range(len(self.bams)):
			roc_fn = self.bams[i].roc_fn()
			t_roc_fn = os.path.basename(roc_fn)

			gp_fn = self.bams[i].gp_fn()
			t_gp_fn = os.path.basename(gp_fn)

			svg_fn = self.bams[i].svg_fn()
			t_svg_fn = os.path.basename(svg_fn)

			add_file_to_tar(tar, roc_fn, t_roc_fn)
			add_file_to_tar(tar, gp_fn, t_gp_fn,
				lambda x: strip_lines(x.replace(roc_fn, t_roc_fn).replace(svg_fn, t_svg_fn)))

		gp_fn = self._gp_fn
		t_gp_fn = os.path.basename(gp_fn)
		svg_dir = os.path.join(self.panel_dir, "graphics") + "/"
		roc_dir = os.path.join(self.panel_dir, "roc") + "/"
		add_file_to_tar(tar, gp_fn, t_gp_fn, lambda x: strip_lines(x.replace(svg_dir, "").replace(roc_dir, "")))

		makefile = [
			".PHONY: all",
			"all:",
			"\tgnuplot *.gp",
			"clean:",
			"\trm -f *.svg",
			"",
		]
		add_text_to_tar(tar, "Makefile", "\n".join(makefile))

	# for svg_fn in self._svg_fns:
	# 	t_roc_fn=os.path.basename(svg_fn)
	# 	add_to_tar(tar,svg_fn,t_svg_fn)

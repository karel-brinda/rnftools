import rnftools
import os
import textwrap

import bs4

from . import DEFAULT_ALLOWED_DELTA
from . import _default_gp_style_func


##############
##############
### REPORT ###
##############
##############
class Report:
	"""Class for an entire report.

	Args:
		name (str): Name of the report (name of output file and dir).
		title (str): Title of the report (if None, then name is used).
		description (str): Description of the report.
		bam_dirs (list of str): Directories with BAM files (this option is mutually exclusive with the 'panels' option).
		panels (list of dicts): Advanced configuration for panels (list of dictionaries with the following keys: 'bam_dir' (mandatory); 'panel_dir', 'name', 'title' (optional)).
		allowed_delta (int): Tolerance of difference of coordinates between true (i.e., expected) alignment and real alignment (very important parameter!).
		default_x_run ((float,float)): Range for x-axis in GnuPlot plots.
		default_y_run ((float,float)): Range for y-axis in GnuPlot plots.
		default_pdf_size_cm ((float,float)): Legacy parameter (does not have any effect).
		default_svg_size_px ((int,int)): Size of SVG picture.
		render_pdf_method (str): Method for svg42pdf to render PDF (None / 'any' (default) / 'cairo' / 'reportlab' / 'inkscape' / 'imagemagick' / 'wkhtmltopdf').
		keep_intermediate_files (bool): Keep files created in intermediate steps during evaluation.
		compress_intermediate_files (bool): Compress files created in intermediate steps during evaluation.
		default_x_axis (str): Values on x-axis, e.g., "({m}+{w})/({M}+{m}+{w})".
		default_x_label (str): Label on x-axis.
		gp_style_func (function(i, nb)): Function assigning GnuPlot styles for overall graphs. Arguments: i: 0-based id of curve, nb: number of curves.		
	"""

	# todo: describe format of panels


	def __init__(
			self,
			name,
			title=None,
			description="",
			allowed_delta=DEFAULT_ALLOWED_DELTA,
			bam_dirs=None,
			panels=None,
			default_x_run=(0.00001, 1.0),
			default_y_run=(60, 100),
			default_pdf_size_cm=(10, 10),
			default_svg_size_px=(640, 640),
			render_pdf_method='any',
			keep_intermediate_files=False,
			compress_intermediate_files=True,
			default_x_axis="({m}+{w})/({M}+{m}+{w})",
			default_x_label="FDR in mapping {{/:Italic(#wrongly mapped reads / #mapped reads)}}  ",
			gp_style_func=_default_gp_style_func,
	):

		self._gp_style_func = gp_style_func

		rnftools.lavender.add_report(self)

		self.name = name
		self.report_dir = self.name
		self.title = self.name if title == None else title
		self.description = description

		self.default_x_run = self._load_x_run(default_x_run)
		self.default_y_run = self._load_y_run(default_y_run)
		self.default_svg_size_px = self._load_svg_size_px(default_svg_size_px)
		self.default_x_label = default_x_label

		self.render_pdf_method = render_pdf_method

		self.allowed_delta = int(allowed_delta)
		assert 0 <= allowed_delta

		self._html_fn = name + ".html"

		assert bam_dirs is None or panels is None, "Panels can be specified using bam_dirs or panels, but not both at the same time."

		self.panels = []

		if bam_dirs is not None:
			assert hasattr(bam_dirs, '__iter__'), "bamdirs should be iterable (list, tuple,  etc.)"

			# assert isinstance(bam_dirs,collections.iterable)
			# assert isinstance(bamdirs, basestring)

			self.panels = [
				rnftools.lavender.Panel(
					bam_dir=bam_dirs[i],
					panel_dir=os.path.join(self.report_dir, str(i)),
					report=self,
					name=str(i),
					title="dir {}".format(i),
					keep_intermediate_files=keep_intermediate_files,
					compress_intermediate_files=compress_intermediate_files,
					default_x_axis=default_x_axis,
					default_x_label=default_x_label,
					gp_style_func=self._gp_style_func,
					render_pdf_method=self.render_pdf_method,
				)
				for i in range(len(bam_dirs))
			]

		if panels is not None:

			for i, panel_dict in enumerate(panels):

				bam_dir = panel_dict["bam_dir"]

				try:
					panel_dir = panel_dict["panel_dir"]
				except KeyError:
					panel_dir = os.path.join(self.report_dir, str(i))

				try:
					panel_name = panel_dict["name"]
				except KeyError:
					panel_name = "panel_{}".format(i)

				try:
					panel_title = panel_dict["title"]
				except KeyError:
					panel_title = "dir {}".format(i)

				self.panels.append(
					rnftools.lavender.Panel(
						bam_dir=bam_dir,
						panel_dir=panel_dir,
						report=self,
						name=panel_name,
						title=panel_title,
						keep_intermediate_files=keep_intermediate_files,
						compress_intermediate_files=compress_intermediate_files,
						default_x_axis=default_x_axis,
						default_x_label=default_x_label,
						gp_style_func=self._gp_style_func,
						render_pdf_method=self.render_pdf_method,
					)
				)

		rnftools.lavender.add_input(self._html_fn)

		# first graph
		self.add_graph(
			"({M}+{m}+{w}+{P}+{x})/{all}",
			title="Mapped reads in all reads",
			y_label="#mapped reads / #reads"
		)
		self.add_graph(
			"{M}/({M}+{m}+{w}+{P}+{x})",
			title="Correctly mapped reads in all mapped reads",
			y_label="#correctly mapped reads / #mapped reads",
		)

		self.add_graph(
			"{M}/({M}+{w}+{x}+{u}+{t}+{P})",
			title="Correctly mapped reads in all reads which should be mapped",
			y_label="#correctly mapped reads / #reads which should be mapped",
		)

		self.add_graph(
			"({U}+{T})/({u}+{U}+{t}+{T})",
			title="Correctly unmapped reads in all unmapped reads",
			y_label="#correctly unmapped reads / #unmapped reads",
		)

		self.add_graph(
			"({U}+{T})/({U}+{T}+{m})",
			title="Correctly unmapped reads in all reads which should be unmapped",
			y_label="#correctly unmapped reads / #reads which should be unmapped",
		)

	def add_graph(
			self,
			y,
			x_label=None,
			y_label="",
			title="",
			x_run=None,
			y_run=None,
			svg_size_px=None,
			key_position="bottom right",
	):
		"""
		Add a new graph to the overlap report.

		Args:
			y (str): Value plotted on y-axis.
			x_label (str): Label on x-axis.
			y_label (str): Label on y-axis.
			title (str): Title of the plot.
			x_run ((float,float)): x-range.
			y_run ((int,int)): y-rang.
			svg_size_px ((int,int): Size of SVG image in pixels.
			key_position (str): GnuPlot position of the legend.
		"""

		if x_run is None:
			x_run = self.default_x_run
		if y_run is None:
			y_run = self.default_y_run
		if svg_size_px is None:
			svg_size_px = self.default_svg_size_px

		for panel in self.panels:
			x_run = self._load_x_run(x_run)
			y_run = self._load_y_run(y_run)
			svg_size_px = self._load_svg_size_px(svg_size_px)
			panel.add_graph(
				y=y,
				x_run=x_run,
				y_run=y_run,
				svg_size_px=svg_size_px,
				y_label=y_label,
				x_label=x_label if x_label is not None else self.default_x_label,
				title=title,
				key_position=key_position,
			)

	def get_report_dir(self):
		"""Get directory report's auxiliary files."""

		return self.report_dir

	def clean(self):
		"""Remove all temporary files."""

		rnftools.utils.shell('rm -fR "{}" "{}"'.format(self.report_dir, self._html_fn))

	def get_panels(self):
		"""Get all contained panels."""

		return self.panels

	######################################
	######################################

	def html_fn(self):
		"""Get name of the HTML file of the report."""

		return self._html_fn

	######################################
	######################################

	def create_html(self):
		"""Create HTML report."""

		html_table = ""
		columns = [panel.get_html_column() for panel in self.panels]
		trs = len(columns[0])
		html_table += os.linesep.join([
			"<tr>{}</tr>".format(
				"".join(
					[
						"<td>{}</td>".format(
							columns[col][row]
						) for col in range(len(columns))
					]
				)
			)
			for row in range(trs)
		])

		with open(self._html_fn, "w+") as f:
			css_src = textwrap.dedent("""\
					.main_table                       {border-collapse:collapse;margin-top:15px;}
					td                                {border: solid #aaaaff 1px;padding:4px;vertical-alignment:top;}
					colgroup, thead                   {border: solid black 2px;padding 2px;}
					.configuration                    {font-size:85%;}
					.configuration, .configuration *  {margin:0;padding:0;}
					.formats                          {text-align:center;margin:20px 0px;}
					img                               {min-width:640px}
			""")

			html_src = """<!DOCTYPE html>
			<html>
			<head>
				<meta charset="UTF-8" />
				<title>{title}</title>
				<style type="text/css">
				{css}
				</style>
			</head>
			<body>
				<h1>{title}</h1>
				<strong>{description}</strong>

				<table class="main_table">
				{html_table}
				</table>

			</body>
			""".format(
				html_table=html_table,
				css=css_src,
				title=self.title,
				description=self.description,
			)

			tidy_html_src = bs4.BeautifulSoup(html_src).prettify()
			f.write(tidy_html_src)

	######################################
	######################################

	@staticmethod
	def _load_x_run(x_run):
		assert len(x_run) == 2
		to_return = [float(x) for x in x_run]
		assert 0.0 < to_return[0] <= 1.0
		return to_return

	@staticmethod
	def _load_y_run(y_run):
		assert len(y_run) == 2
		to_return = [float(x) for x in y_run]
		assert 0.0 <= to_return[0] <= 100.0
		assert 0.0 <= to_return[1] <= 100.0
		return to_return

	@staticmethod
	def _load_svg_size_px(svg_size_px):
		assert len(svg_size_px) == 2
		to_return = [int(x) for x in svg_size_px]
		assert 0 <= to_return[0]
		assert 0 <= to_return[1]
		return to_return

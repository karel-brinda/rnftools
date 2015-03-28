import rnftools
import snakemake
import os


##############
##############
### REPORT ###
##############
##############
class Report:
	"""Class for an entire report."""

	def __init__(
		self,
		name,
		bam_dirs=[],
		allowed_delta=5,
		default_plot_x_run=(0.00001,1.0),
		default_plot_y_run=(60,100),
		default_plot_pdf_size_cm=(10,10),
		default_plot_svg_size_px=(640,640),
		keep_aci=False,
		compress_aci=True,
		default_x_axis="({m}+{w})/({M}+{m}+{w})",
		):
		"""

		:param name: Name of the report.
		:type  name: str
		:param bam_dirs: Directories with BAM files.
		:type  bam_dirs: list of str
		:param allowed_delta: Tolerance of difference between positions in assessing correct alignments (very important parameter!).
		:type  allowed_delta: int
		:param default_plot_x_run: Range for x-axis in GnuPlot plots.
		:type  default_plot_x_run: (float,float)
		:param default_plot_y_run: Range for y-axis in GnuPlot plots.
		:type  default_plot_y_run: (float,float)
		:param default_plot_pdf_size_cm: Size of PDF page.
		:type  default_plot_pdf_size_cm: (float,float)
		:param default_plot_svg_size_px: Size of SVG picture.
		:type  default_plot_svg_size_px: (int,int)

		"""

		rnftools.lavender.add_report(self)

		self.name = name
		self.report_dir = self.name

		self.default_plot_x_run=self._load_plot_x_run(default_plot_x_run)
		self.default_plot_y_run=self._load_plot_y_run(default_plot_y_run)
		self.default_plot_pdf_size_cm=self._load_pdf_size_cm(default_plot_pdf_size_cm)
		self.default_plot_svg_size_px=self._load_svg_size_px(default_plot_svg_size_px)

		self.allowed_delta=int(allowed_delta)
		assert 0 <= allowed_delta

		self._html_fn = name+".html"
		self.panels = [
				rnftools.lavender.Panel(
					bam_dir=bam_dirs[i],
					panel_dir=os.path.join(self.report_dir,str(i)),
					report=self,
					name=str(i),
					keep_aci=keep_aci,
					compress_aci=compress_aci,
					default_x_axis=default_x_axis,
				)
				for i in range(len(bam_dirs))
			]

		rnftools.lavender.add_input(self._html_fn)

		# first graph
		self.add_graph("({M}+{m}+{w})/{all}")

	def add_graph(self,
				y,
				ylabel="",
				plot_x_run=None,
				plot_y_run=None,
				plot_pdf_size_cm=None,
				plot_svg_size_px=None,
			):

		if plot_x_run==None:
			plot_x_run=self.default_plot_x_run
		if plot_y_run==None:
			plot_y_run=self.default_plot_y_run
		if plot_pdf_size_cm==None:
			plot_pdf_size_cm=self.default_plot_pdf_size_cm
		if plot_svg_size_px==None:
			plot_svg_size_px=self.default_plot_svg_size_px

		for panel in self.panels:
			plot_x_run=self._load_plot_x_run(plot_x_run)
			plot_y_run=self._load_plot_y_run(plot_y_run)
			plot_pdf_size_cm=self._load_pdf_size_cm(plot_pdf_size_cm)
			plot_svg_size_px=self._load_svg_size_px(plot_svg_size_px)
			panel.add_graph(
					y=y,
					plot_x_run=plot_x_run,
					plot_y_run=plot_y_run,
					plot_pdf_size_cm=plot_pdf_size_cm,
					plot_svg_size_px=plot_svg_size_px,
					ylabel=ylabel,
				)

	def get_report_dir(self):
		"""Get directory report's auxiliary files."""
		
		return self.report_dir

	def clean(self):
		"""Remove all temporary files."""

		snakemake.shell('rm -fR "{}" "{}"'.format(self.report_dir,self._html_fn))

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
		html_table+=os.linesep.join([
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

		with open(self._html_fn,"w+") as f:
			html_src="""<!DOCTYPE html>
			<html>
			<head>
				<meta charset="UTF-8" />
				<title>{name}</title>
				<style type="text/css">
					table                             {{border-collapse:collapse}}
					td                                {{border: solid #aaaaff 1px;padding:4px;vertical-alignment:top}}
					colgroup, thead                   {{border: solid black 2px;padding 2px}}
					td img                            {{width:100%}}
					.configuration                    {{font-size:85%}}
					.configuration, .configuration *  {{margin:0;padding:0}}
				</style>
			</head>
			<body>
				<h1>{name}</h1>
				<table>
				{html_table}
				</table>
			</body>
			""".format(
					html_table=html_table,
					name=self.name
				)
			f.write(html_src)

	######################################
	######################################

	@staticmethod
	def _load_plot_x_run(plot_x_run):
		assert len(plot_x_run)==2
		to_return=[float(x) for x in plot_x_run]
		assert 0< to_return[0] and to_return[0]<=1.0
		return to_return

	@staticmethod
	def _load_plot_y_run(plot_y_run):
		assert len(plot_y_run)==2
		to_return=[float(x) for x in plot_y_run]
		assert 0<=to_return[0] and to_return[0]<=100
		assert 0<=to_return[1] and to_return[1]<=100
		return to_return

	@staticmethod
	def _load_pdf_size_cm(pdf_size_cm):
		assert len(pdf_size_cm)==2
		to_return=[float(x) for x in pdf_size_cm]
		assert 0<=to_return[0]
		assert 0<=to_return[1]
		return to_return

	@staticmethod
	def _load_svg_size_px(svg_size_px):
		assert len(svg_size_px)==2
		to_return=[int(x) for x in svg_size_px]
		assert 0<=to_return[0]
		assert 0<=to_return[1]
		return to_return

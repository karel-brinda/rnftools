import rnftools
import snakemake
import os


__all__=["Report"]

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
		plot_x_run=(0.00001,1.0),
		plot_y_run=(60,100),
		plot_pdf_size_cm=(10,10),
		plot_svg_size=(640,640),
		):
		"""

		:param name: Name of the report.
		:type  name: str
		:param bam_dirs: Directories with BAM files.
		:type  bam_dirs: list of str
		:param allowed_delta: Tolerance of difference between positions in assessing correct alignments.
		:type  allowed_delta: int
		:param plot_x_run: Range for x-axes in GnuPlot plots.
		:type  plot_x_run: (float,float)
		:param plot_y_run: Range for y-axes in GnuPlot plots.
		:type  plot_y_run: (float,float)
		:param plot_pdf_size_cm: Size of PDF page.
		:type  plot_pdf_size_cm: (float,float)
		:param plot_svg_size: Size of SVG picture.
		:type  plot_svg_size: (int,int)

		"""

		rnftools.lavender.add_report(self)

		self.name = name
		self.report_dir = self.name

		self.plot_x_run=[float(x) for x in plot_x_run]
		self.plot_y_run=[float(x) for x in plot_y_run]
		self.plot_pdf_size_cm=[float(x) for x in plot_pdf_size_cm]
		self.plot_svg_size=[int(x) for x in plot_svg_size]

		assert 0< self.plot_x_run[0] and self.plot_x_run[0]<=1.0
		assert 0< self.plot_x_run[1] and self.plot_x_run[1]<=1.0
		assert 0<=self.plot_y_run[0] and self.plot_y_run[0]<=100
		assert 0<=self.plot_y_run[1] and self.plot_y_run[1]<=100
		assert 0<=self.plot_pdf_size_cm[0] 
		assert 0<=self.plot_pdf_size_cm[1]  
		assert 0<=self.plot_svg_size[0]
		assert 0<=self.plot_svg_size[1]

		self.allowed_delta=int(allowed_delta)
		assert 0 <= allowed_delta


		self._html_fn = name+".html"
		self.panels = [
				rnftools.lavender.Panel(
					bam_dir=bam_dirs[i],
					panel_dir=os.path.join(self.report_dir,str(i)),
					report=self,
					name=str(i)
				)
				for i in range(len(bam_dirs))
			]

		rnftools.lavender.add_input(self._html_fn)

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


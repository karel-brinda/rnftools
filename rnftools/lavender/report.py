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

	def __init__(self,name,bam_dirs=[]):
		"""

		:param name: Name of the report.
		:type  name: str
		:param bam_dirs: Directories with BAM files.
		:type  bam_dirs: list of str

		"""

		rnftools.lavender.add_report(self)

		self.name = name
		self.report_dir = self.name

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


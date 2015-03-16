# -*- coding: utf-8 -*-

import rnftools.lavender

import smbl.prog
import snakemake
import os
import sys
import pysam


__all__ = ["Bam"]

MAXIMAL_MAPPING_QUALITY=255

###########
###########
### BAM ###
###########
###########
class Bam:
	"""Class for a single BAM file."""

	def __init__(self,panel,bam_fn,name):
		"""

		:param panel: Panel containing this BAM file.
		:type  panel: rnftools.lavender.Panel
		:param bam_fn: BAM filename.
		:type  bam_fn: str
		:param name: Name for this report.
		:type  name: str

		"""

		self.panel=panel
		self.report=panel.get_report()
		self.name=name

		self._bam_fn  = bam_fn
		self._gp_fn   = os.path.join(self.panel.get_panel_dir(),"gp",self.name+".gp")
		self._html_fn = os.path.join(self.panel.get_panel_dir(),"html",self.name+".html")
		self._aci_fn  = os.path.join(self.panel.get_panel_dir(),"aci",self.name+".aci")
		self._roc_fn  = os.path.join(self.panel.get_panel_dir(),"roc",self.name+".roc")
		self._svg_fn  = os.path.join(self.panel.get_panel_dir(),"svg",self.name+".svg")
		self._pdf_fn  = os.path.join(self.panel.get_panel_dir(),"pdf",self.name+".pdf")

		self.bam_id=len(rnftools.lavender.bams())
		rnftools.lavender.add_bam(self)

	def get_name(self):
		"""Get name associated with the BAM."""
		return self.name

	def get_required_files(self):
		"""Get names of all files required to complete the report."""
		return [self._svg_fn, self._html_fn]

	######################################
	######################################

	def bam_fn(self):
		"""Get name of the BAM file."""
		return self._bam_fn

	def gp_fn(self):
		"""Get name of the GP file."""
		return self._gp_fn

	def html_fn(self):
		"""Get name of the HTML report."""
		return self._html_fn

	def aci_fn(self):
		"""Get name of the ACIfile."""
		return self._aci_fn

	def roc_fn(self):
		"""Get name of the ROC file."""
		return self._roc_fn

	def svg_fn(self):
		"""Get name of the SVG file."""
		return self._svg_fn

	def pdf_fn(self):
		"""Get name of the PDF file."""
		return self._pdf_fn

	######################################
	######################################

	def create_gp(self):
		"""Create a GnuPlot file for this BAM file."""

		with open(self._gp_fn,"w+") as f:
			gp_content="""
				set key spacing 0.8 opaque width -3

				set x2lab "false positive rate\\n(#wrong mappings / #mapped)"
				set log x
				set log x2

				set format x "10^{{%L}}"
				set format x2 "10^{{%L}}"
				set xran  [{xran}]
				set x2ran [{xran}]
				set x2tics


				set ylab "part of all reads (%)"

				set format y "%g %%"
				set yran [{yran}]

				set pointsize 1.5

				set grid ytics lc rgb "#777777" lw 1 lt 0 front
				set grid xtics lc rgb "#777777" lw 1 lt 0 front

				set datafile separator "\t"
				set palette negative

				set termin svg size {svg_size} enhanced
				set out "{svg_fn}"
				{plot}


				set termin pdf enhanced size {pdf_size} enhanced font 'Arial,12'
				set out "{pdf_fn}"
				{plot}

			""".format(
				svg_fn=self._svg_fn,
				pdf_fn=self._pdf_fn,
				xran="{:.10f}:{:.10f}".format(self.report.plot_x_run[0],self.report.plot_x_run[1]),
				yran="{:.10f}:{:.10f}".format(self.report.plot_y_run[0],self.report.plot_y_run[1]),
				svg_size="{},{}".format(self.report.plot_svg_size[0],self.report.plot_svg_size[1]),
				pdf_size="{:.10f}cm,{:.10f}cm".format(self.report.plot_pdf_size_cm[0],self.report.plot_pdf_size_cm[1]),
				plot="""
					plot\\
						"{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2+$3+$4+$8+$7+$6+$5)*100/$10) lt rgb "violet" with filledcurve x1 title 'Unmapped correctly',\\
						"{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2+$3+$4+$8+$7+$6)*100/$10) lt rgb "red" with filledcurve x1 title 'Unmapped incorrectly',\\
						"{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2+$3+$4+$8+$7)*100/$10) lt rgb "green" with filledcurve x1 title 'Thresholded',\\
						"{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2+$3+$4+$8)*100/$10) lt rgb "yellow" with filledcurve x1 title 'Multimapped',\\
						"{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2+$3+$4)*100/$10) lt rgb "gray" with filledcurve x1 title 'Mapped, should be unmapped',\\
						"{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2+$3)*100/$10) lt rgb "black" with filledcurve x1 title 'Mapped to wrong position',\\
						"{roc_fn}" using (( ($3+$4) / ($2+$3+$4) )):(($2)*100/$10) lt rgb "blue" with filledcurve x1 title 'Mapped correctly',\\

				""".format(roc_fn=self._roc_fn)
			)
			f.write(gp_content)

	def create_html(self):
		"""Create a HTML page for this BAM file."""

		roc_dicts = []
		with open(self._roc_fn,"r") as f:
			for line in f:
				line=line.strip()
				if line!="" and line[0]!="#":
					(q,M,w,m,U,u,t,p,x,a)=line.split("\t")
					roc_dict = {
						"q":int(q),
						"M":int(M),
						"w":int(w),
						"m":int(m),
						"U":int(U),
						"u":int(u),
						"t":int(t),
						"+":int(p),
						"x":int(x),
						"a":int(a)
					}
					roc_dicts.append(roc_dict)
		tbody = os.linesep.join([
				"""
					<tr>
						<td>        {quality}                   </td>
						<td>        {mapped}                    </td>
						<td><small> {mapped_proc:.2f}   </small></td>
						<td>        {M}                         </td>
						<td><small> {M_proc:.2f}        </small></td>
						<td>        {w}                         </td>
						<td><small> {w_proc:.2f}        </small></td>
						<td>        {m}                         </td>
						<td><small> {m_proc:.2f}        </small></td>
						<td>        {unmapped}                  </td>
						<td><small> {unmapped_proc:.2f} </small></td>
						<td>        {U}                         </td>
						<td><small> {U_proc:.2f}        </small></td>
						<td>        {u}                         </td>
						<td><small> {u_proc:.2f}        </small></td>
						<td>        {unused}                    </td>
						<td><small> {unused_proc:.2f}   </small></td>
						<td>        {t}                         </td>
						<td><small> {t_proc:.2f}        </small></td>
						<td>        {p}                         </td>
						<td><small> {p_proc:.2f}        </small></td>
						<td>        {x}                         </td>
						<td><small> {x_proc:.2f}        </small></td>
						<td>        {sum}                       </td>
						<td>        {prec_proc:.3f}              </td>
						<td>        {sens_proc:.3f}             </td>
					</tr>
				""".format(
						quality       = roc_dict["q"],
						mapped        = roc_dict["M"]+roc_dict["w"]+roc_dict["m"],
						mapped_proc   = 100.0*(roc_dict["M"]+roc_dict["w"]+roc_dict["m"])/roc_dict["a"],
						M             = roc_dict["M"],
						M_proc        = 100.0*(roc_dict["M"])/roc_dict["a"],
						w             = roc_dict["w"],
						w_proc        = 100.0*(roc_dict["w"])/roc_dict["a"],
						m             = roc_dict["m"],
						m_proc        = 100.0*(roc_dict["m"])/roc_dict["a"],
						unmapped      = roc_dict["U"]+roc_dict["u"],
						unmapped_proc = 100.0*(roc_dict["U"]+roc_dict["u"])/roc_dict["a"],
						U             = roc_dict["U"],
						U_proc        = 100.0*(roc_dict["U"])/roc_dict["a"],
						u             = roc_dict["u"],
						u_proc        = 100.0*(roc_dict["u"])/roc_dict["a"],
						unused        = roc_dict["t"]+roc_dict["+"]+roc_dict["x"],
						unused_proc   = 100.0*(roc_dict["t"]+roc_dict["+"]+roc_dict["x"])/roc_dict["a"],
						t             = roc_dict["t"],
						t_proc        = 100.0*(roc_dict["t"])/roc_dict["a"],
						p             = roc_dict["+"],
						p_proc        = 100.0*(roc_dict["+"])/roc_dict["a"],
						x             = roc_dict["x"],
						x_proc        = 100.0*(roc_dict["x"])/roc_dict["a"],
						sum           = roc_dict["a"],
						sens_proc     = 100.0*(roc_dict["M"]+roc_dict["w"]+roc_dict["m"])/roc_dict["a"],
						prec_proc      = 100.0*(roc_dict["M"])/(roc_dict["M"]+roc_dict["w"]+roc_dict["m"]) if (roc_dict["M"]+roc_dict["w"]+roc_dict["m"]) != 0 else 0,
					)

				for roc_dict in roc_dicts
			])


		with open(self._html_fn,"w+") as f:
			html_src="""<!DOCTYPE html>
			<html>
			<head>
				<meta charset="UTF-8" />
				<title>{name}</title>
				<style type="text/css">
					table            {{border-collapse:collapse}}
					td               {{border: solid #aaaaff 1px;padding:4px;text-align:right}}
					colgroup, thead  {{border: solid black 2px;padding 2px}}
					.link_to_top     {{font-size:10pt}}
				</style>
			</head>
			<body>
				<h1 id="top">{name}</h1>

				<p>
					<a href="#roctable">ROC table</a> -
					<a href="#graphs">Graphs</a>
				</p>

				<h2 id="roctable">
					ROC table
					<span class="link_to_top">
						[<a href="#top">Top of the page</a>]
						[<a href="{homepage}">Main report</a>]
					</span>
				</h2>

				<p>
					<strong>M</strong>: mapped correctly,
					<strong>w</strong>: mapped to wrong position,
					<strong>m</strong>: mapped but should be unmapped,
					<strong>U</strong>: unmapped correctly,
					<strong>u</strong>: unmapped but should be mapped,
					<strong>t</strong>: thresholded,
					<strong>+</strong>: multimapped,
					<strong>x</strong>: unknown (read is probably not present in SAM)
				</p>
				<table>
					<colgroup span="1" style="">
					<colgroup span="2" style="background-color:#ddd">
					<colgroup span="6" style="">
					<colgroup span="2" style="background-color:#ddd">
					<colgroup span="4" style="">
					<colgroup span="2" style="background-color:#ddd">
					<colgroup span="6" style="">
					<colgroup span="1" style="background-color:#ddd">
					<colgroup span="2" style="">
					<thead style="font-weight:bold;background-color:#ddddff">
						<tr style="font-weight:bold;background-color:#ddddff">
							<td>q</td>
							<td>mapped</td>
							<td>%</td>
							<td>M</td>
							<td>%</td>
							<td>w</td>
							<td>%</td>
							<td>m</td>
							<td>%</td>
							<td>unmapped</td>
							<td>%</td>
							<td>U</td>
							<td>%</td>
							<td>u</td>
							<td>%</td>
							<td>unused</td>
							<td>%</td>
							<td>t</td>
							<td>%</td>
							<td>+</td>
							<td>%</td>
							<td>x</td>
							<td>%</td>
							<td>sum</td>
							<td>prec. (%)</td>
							<td>sens. (%)</td>
						</tr>
					</thead>
					<tbody>
						{tbody}
					</tbody>
				</table>

				<h2 id="graphs">
					Graphs
					<span class="link_to_top">
						[<a href="#top">Top of this page</a>]
						[<a href="{homepage}">Back to Main report</a>]
					</span>
				</h2>
				<img src="{svg_fn}" />
			</body>
			</html>			
			""".format(
					name=self.name,
					tbody=tbody,
					svg_fn=os.path.relpath(
						self._svg_fn,
						os.path.dirname(self._html_fn)
					),
					homepage=os.path.relpath(
						self.report.html_fn(),
						os.path.dirname(self._html_fn)
					),
				)

			f.write(html_src)

	def create_aci(self):
		"""Create an ACI (intermediate) file for this BAM file.
		This is the function which asses if an alignment is correct
		"""

		diff_thr=self.report.allowed_delta

		with open(self._aci_fn,"w+") as f:
			with pysam.AlignmentFile(self._bam_fn, "rb") as samfile:
				references_dict = {}

				for i in range(len(samfile.references)):
					references_dict[ samfile.references[i] ] = i+1

				f.write("# read name"+os.linesep)
				f.write("# is mapped with quality"+os.linesep)
				f.write("# chr id"+os.linesep)
				f.write("# direction"+os.linesep)
				f.write("# the most left nucleotide"+os.linesep)
				f.write("# the most right nucleotide"+os.linesep)
				f.write("# category of alignment assigned by LAVEnder"+os.linesep)
				f.write("#      M_i    i-th bsegment is correctly mapped"+os.linesep)
				f.write("#      m      segment should be unmapped but it is mapped"+os.linesep)
				f.write("#      w      segment is mapped to a wrong location"+os.linesep)
				f.write("#      U      segment is unmapped and should be unmapped"+os.linesep)
				f.write("#      u      segment is unmapped and should be mapped"+os.linesep)
				f.write("# number of segments"+os.linesep)

				for read in samfile:
					rnf_read = rnftools.rnfformat.Read()
					rnf_read.destringize(read.query_name)

					left = read.reference_start+1
					right = read.reference_end
					chrom_id=references_dict[ samfile.references[read.reference_id] ]

					nb_of_segments=len(rnf_read.segments)

					if rnf_read.segments[0].source==1:
						should_be_mapped=True
					else:
						should_be_mapped=False

					# read is unmapped
					if read.is_unmapped:
						# read should be mapped
						if should_be_mapped:
							category="u"			
						# read should be unmapped
						else:
							category="U"
					# read is mapped
					else:
						# read should be mapped
						if should_be_mapped:
							exists_corresponding_segment=False

							for j in range(len(rnf_read.segments)):
								segment=rnf_read.segments[j]
								if (
									(segment.left==0 or abs(segment.left - left) < diff_thr) and
									(segment.right==0 or abs(segment.right - right) < diff_thr) and
									(segment.left!=0 or segment.right==0) and
									(chrom_id==0 or chrom_id==segment.chr)
								):
									exists_corresponding_segment=True
									segment=str(j+1)
									break

							# read was mapped to correct location
							if exists_corresponding_segment: # exists ok location?
								category="M_"+segment
							# read was mapped to incorrect location
							else:
								category="w"
						# read should be unmapped
						else:
							category="m"


					f.write(
						"\t".join(
							map(str,[
								# read name
								read.query_name,
								# aligned?
								"unmapped" if read.is_unmapped else "mapped_"+str(read.mapping_quality),
								# reference id
								chrom_id,
								# direction
								"R" if read.is_reverse else "F",
								# left
								left,
								# right
								right,
								# assigned category
								category,
								# count of segments
								nb_of_segments
							])
						) + os.linesep
					)

	@staticmethod
	def _vector_of_categories(srs,read_name,parts):
		"""Create vector of categories (voc[q] ... assigned category for given quality level)

		srs ... single read statistics ... for every q ... dictionary
		read_name ... read name
		parts ... number of segments
		"""

		# default value
		vec = ["x" for i in range(MAXIMAL_MAPPING_QUALITY+1)]
		val_err = lambda q : """\
Invalid alignment for read '{}'. Debug info:
	{dbg}
Please contact the author on karel.brinda@gmail.com.
""".format(read_name,dbg=str( [
		srs[q]["M"],
		srs[q]["m"],
		srs[q]["w"],
		srs[q]["U"],
		srs[q]["u"]
	]))
		assert len(srs)<=MAXIMAL_MAPPING_QUALITY+1

		for q in range(len(srs)):

			#assert (
			#		len(set(srs[q]["M"]))+
			#		srs[q]["m"]+
			#		srs[q]["w"]+
			#		srs[q]["U"]+
			#		srs[q]["u"]) <= parts+1, "wrong single read statistics parts={}, q={}, '{}'".format(parts,q,str(srs[q]))


			#####
			# M # - all parts correctly aligned
			#####
			if (
				len(srs[q]["M"])==parts and
				srs[q]["w"]==0 and
				srs[q]["m"]==0 and
				srs[q]["U"]==0 and
				srs[q]["u"]==0
			):
				assert vec[q]=="x",str((q,srs[q]))
				vec[q]="M"

			#####
			# w # - at least one segment is incorrectly aligned
			#####
			if (
				srs[q]["w"]>0
			):
				assert vec[q]=="x",str((q,srs[q]))
				vec[q]="w"

			#####
			# m # - at least one segment was aligned but should not be aligned
			#####
			if(
				srs[q]["w"]==0 and
				srs[q]["m"]>0
			):
				assert vec[q]=="x",str((q,srs[q]))
				vec[q]="m"

			#####
			# U # - all segments should be unaligned but are unaligned
			#####
			if (
				srs[q]["U"]>0 and
				srs[q]["u"]==0 and
				srs[q]["m"]==0 and
				srs[q]["w"]==0 and
				len(srs[q]["M"])==0
			):
				assert vec[q]=="x",str((q,srs[q]))
				vec[q]="U"

			#####
			# u # - at least one segment was unaligned but should be aligned
			#####
			if (
				srs[q]["w"]==0 and
				srs[q]["u"]>0
			):
				assert vec[q]=="x",str((q,srs[q]))
				vec[q]="u"

			#####
			# t # - at least one segment was thresholded
			#####
			if (
				len(srs[q]["M"])!=parts and
				srs[q]["w"]==0 and
				srs[q]["m"]==0 and
				srs[q]["U"]==0 and
				srs[q]["u"]==0
			):
				assert vec[q]=="x",str((q,srs[q]))
				vec[q]="t"

			#####
			# + # - multimapped, M + w + m > parts
			#####

			# this only can remove some older
			if (
				#len(srs[q]["M"])+srs[q]["w"]+srs[q]["m"]>parts+srs[q]["m"]==0 and
				len(srs[q]["M"])+srs[q]["w"]+srs[q]["m"]>parts and
				srs[q]["U"]==0 and
				srs[q]["u"]==0
			):
				#assert vec[q]=="x",str((q,srs[q]))
				vec[q]="+"
				#print("multi")

		return vec


	def create_roc(self):
		"""Create a ROC file for this BAM file.

		:raises: ValueError

		"""

		stats_dicts = [
			{
				"q":i,
				"M":0,
				"w":0,
				"m":0,
				"U":0,
				"u":0,
				"t":0,
				"+":0,
				"x":0
			}
			for i in range(MAXIMAL_MAPPING_QUALITY+1)
		]

		with open(self._aci_fn, "r") as g:

			last_rname=""
			for line in g:
				line=line.strip()
				if line=="" or line[0]=="#":
					continue
				else:
					(rname,mapped,ref,direction,left,right,category,nb_of_segments)=line.split("\t")
					nb_of_segments=int(nb_of_segments)

					#print(rname,last_rname,mapped)
					# new read
					if rname!=last_rname:
						# update
						if last_rname!="":
							voc = self._vector_of_categories(single_reads_statistics,rname,nb_of_segments)
							for q in range(len(voc)):
								stats_dicts[q][voc[q]]+=1
						# nulling
						single_reads_statistics= [
									{
										"U":0,
										"u":0,
										"M":[],
										"m":0,
										"w":0,
									} for i in range(MAXIMAL_MAPPING_QUALITY+1)
								]
						last_rname=rname

					# processing of a segment
					if category=="U":
						for q in range(len(single_reads_statistics)):
							single_reads_statistics[q]["U"]+=1
					elif category=="u":
						for q in range(len(single_reads_statistics)):
							single_reads_statistics[q]["u"]+=1
					else:
						mapping_quality=int(mapped.replace("mapped_",""))
						assert mapping_quality<=MAXIMAL_MAPPING_QUALITY
						if category=="m":
							for q in range(mapping_quality+1):
								single_reads_statistics[q]["m"]+=1
						elif category=="w":
							for q in range(mapping_quality+1):
								single_reads_statistics[q]["w"]+=1
						else:
							assert category[0]=="M", category
							segment_id=int(category.replace("M_",""))
							for q in range(mapping_quality+1):
								single_reads_statistics[q]["M"].append(segment_id)

			# last read
			voc = self._vector_of_categories(single_reads_statistics,rname,nb_of_segments)
			for q in range(len(voc)):
				stats_dicts[q][voc[q]]+=1


			with open(self._roc_fn, "w+") as f:
				f.write("# Numbers of reads in several categories in dependence"+os.linesep)
				f.write("# on the applied threshold on mapping quality q"+os.linesep)
				f.write("# "+os.linesep)
				f.write("# Categories:"+os.linesep)
				f.write("#        M: Mapped correctly."+os.linesep)
				f.write("#        w: Mapped to a wrong position."+os.linesep)
				f.write("#        m: Mapped but should be unmapped."+os.linesep)
				f.write("#        U: Unmapped and should be unmapped."+os.linesep)
				f.write("#        u: Unmapped but should be mapped."+os.linesep)
				f.write("#        t: Thresholded."+os.linesep)
				f.write("#        +: Multimapped."+os.linesep)
				f.write("#        x: Unknown."+os.linesep)
				f.write("#"+os.linesep)
				f.write("# q\tM\tw\tm\tU\tu\tt\t+\tx\tall"+os.linesep)

				l_numbers = []
				for line in stats_dicts:
					numbers = [line["M"],line["w"],line["m"],line["U"],line["u"],line["t"],line["+"],line["x"]]
					if numbers != l_numbers:
						f.write("\t".join(
								[str(line["q"])] + list(map(str,numbers)) + [str(sum(numbers))]
							)+os.linesep)
					l_numbers=numbers


	def create_graphics(self):
		"""Create images related to this BAM file using GnuPlot."""

		snakemake.shell('"{}" "{}"'.format(smbl.prog.GNUPLOT5,self._gp_fn))

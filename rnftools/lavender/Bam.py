# -*- coding: utf-8 -*-

import rnftools.lavender

import smbl.prog
import snakemake
import os
import sys
import pysam
import gzip

###########
###########
### BAM ###
###########
###########
class Bam:
	"""Class for a BAM file.

	Args:
		panel (rnftools.lavender.Panel): Panel containing this BAM file.
		bam_fn (str): BAM filename.
		name (str): Name for this report.
		keep_intermediate_files (bool): Keep files created in intermediate steps during evaluation.
		compress_intermediate_files (bool): Compress files created in intermediate steps during evaluation.
		default_x_axis (str): Values on x-axis, e.g., "({m}+{w})/({M}+{m}+{w})".
		default_x_label (str): Label on x-axis.

	"""

	def __init__(self,
				panel,
				bam_fn,
				name,
				keep_intermediate_files,
				compress_intermediate_files,
				default_x_axis,
				default_x_label,
			):

		self.panel=panel
		self.report=panel.get_report()
		self.name=name

		self.keep_intermediate_files=keep_intermediate_files
		self.compress_intermediate_files=compress_intermediate_files
		self.default_x_axis=default_x_axis
		self.default_x_label=default_x_label

		self._bam_fn  = bam_fn
		self._gp_fn   = os.path.join(self.panel.get_panel_dir(),"gp",self.name+".gp")
		self._html_fn = os.path.join(self.panel.get_panel_dir(),"html",self.name+".html")
		if compress_intermediate_files:
			self._es_fn  = os.path.join(self.panel.get_panel_dir(),"es",self.name+".es.gz")
			self._et_fn  = os.path.join(self.panel.get_panel_dir(),"et",self.name+".et.gz")
		else:
			self._es_fn  = os.path.join(self.panel.get_panel_dir(),"es",self.name+".es")
			self._et_fn  = os.path.join(self.panel.get_panel_dir(),"et",self.name+".et")
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

	def es_fn(self):
		"""Get name of the es file."""
		return self._es_fn

	def et_fn(self):
		"""Get name of the et file."""
		return self._et_fn

	def roc_fn(self):
		"""Get name of the ROC file."""
		return self._roc_fn

	def svg_fn(self):
		"""Get name of the SVG file."""
		return self._svg_fn

	def pdf_fn(self):
		"""Get name of the PDF file."""
		return self._pdf_fn


	############################
	############################
	##
	## ES
	##
	############################
	############################

	@staticmethod
	def bam2es(
				bam_fn,
				es_fo,
				allowed_delta,
			):
		"""Convert BAM file to ES file.

		Args:
			bam_fn (str): File name of the BAM file.
			bam_fo (file): File object of the ES file.
			allowed_delta (int): Maximal allowed coordinates difference for correct reads.
		"""

		es_fo.write("# RN:   read name"+os.linesep)
		es_fo.write("# Q:    is mapped with quality"+os.linesep)
		es_fo.write("# Chr:  chr id"+os.linesep)
		es_fo.write("# D:    direction"+os.linesep)
		es_fo.write("# L:    leftmost nucleotide"+os.linesep)
		es_fo.write("# R:    rightmost nucleotide"+os.linesep)
		es_fo.write("# Cat:  category of alignment assigned by LAVEnder"+os.linesep)
		es_fo.write("#         M_i    i-th segment is correctly mapped"+os.linesep)
		es_fo.write("#         m      segment should be unmapped but it is mapped"+os.linesep)
		es_fo.write("#         w      segment is mapped to a wrong location"+os.linesep)
		es_fo.write("#         U      segment is unmapped and should be unmapped"+os.linesep)
		es_fo.write("#         u      segment is unmapped and should be mapped"+os.linesep)
		es_fo.write("# Segs: number of segments"+os.linesep)
		es_fo.write("# "+os.linesep)
		es_fo.write("# RN\tQ\tChr\tD\tL\tR\tCat\tSegs"+os.linesep)

		with pysam.AlignmentFile(bam_fn, "rb") as sam:
			references_dict = {}

			for i in range(len(sam.references)):
				references_dict[ sam.references[i] ] = i+1

			for read in sam:
				rnf_read_tuple = rnftools.rnfformat.ReadTuple()
				rnf_read_tuple.destringize(read.query_name)

				left = read.reference_start+1
				right = read.reference_end
				chrom_id=references_dict[ sam.references[read.reference_id] ]

				nb_of_segments=len(rnf_read_tuple.segments)

				if rnf_read_tuple.segments[0].genome_id==1:
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

						for j in range(len(rnf_read_tuple.segments)):
							segment=rnf_read_tuple.segments[j]
							if (
								(segment.left==0 or abs(segment.left - left) <= allowed_delta) and
								(segment.right==0 or abs(segment.right - right) <= allowed_delta) and
								(segment.left!=0 or segment.right==0) and
								(chrom_id==0 or chrom_id==segment.chr_id)
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

				es_fo.write(
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

	def create_es(self):
		"""Create an ES (intermediate) file for this BAM file.
		This is the function which asses if an alignment is correct
		"""

		with (gzip.open(self._es_fn,"tw+") if self.compress_intermediate_files else open(self._es_fn,"w+")) as es_fo:
			self.bam2es(
					bam_fn=self._bam_fn,
					es_fo=es_fo,
					allowed_delta=self.report.allowed_delta,
				)


	############################
	############################
	##
	## ET
	##
	############################
	############################


	@staticmethod
	def _vector_of_categories(srs,read_tuple_name,parts):
		"""Create vector of categories (voc[q] ... assigned category for given quality level)

		srs ... single read statistics ... for every q ... dictionary
		read_tuple_name ... read name
		parts ... number of segments
		"""

		# default value
		vec = ["x" for i in range(rnftools.lavender.MAXIMAL_MAPPING_QUALITY+1)]
		assert len(srs)<=rnftools.lavender.MAXIMAL_MAPPING_QUALITY+1,srs

		should_be_mapped=bool(srs[0]["m"]+srs[0]["U"]==0)

		for q in range(len(srs)):

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
				srs[q]["u"]==0 and
				srs[q]["t"]>0
			):
				assert vec[q]=="x",str((q,srs[q]))
				vec[q]="t"

			#####
			# T # - at least one segment was thresholded
			#####
			if (
				len(srs[q]["M"])!=parts and
				srs[q]["w"]==0 and
				srs[q]["m"]==0 and
				srs[q]["U"]==0 and
				srs[q]["u"]==0 and
				srs[q]["T"]>0
			):
				assert vec[q]=="x",str((q,srs[q]))
				vec[q]="T"

			#####
			# P # - multimapped, M + w + m > parts
			#####

			# this only can rewrite some older
			if (
				#len(srs[q]["M"])+srs[q]["w"]+srs[q]["m"]>parts+srs[q]["m"]==0 and
				len(srs[q]["M"])+srs[q]["w"]+srs[q]["m"]>parts and
				srs[q]["U"]==0 and
				srs[q]["u"]==0
			):
				#assert vec[q]=="x",str((q,srs[q]))
				vec[q]="P"


			######
			## x # - unrecognized - print details
			######
			#if vec[q]=="x":
			#	smbl.messages.message(
			#		" ".join(
			#			[
			#				"Unrecognized category for alignment of read '{}'.".format(read_tuple_name),
			#				"Quality level: {}.".format(q),
			#				"Debug info: '{}'.".format(str(
			#						[
			#							srs[q]["M"],
			#							srs[q]["m"],
			#							srs[q]["w"],
			#							srs[q]["U"],
			#							srs[q]["u"],
			#							srs[q]["T"],
			#							srs[q]["t"],
			#						]
			#					)),
			#			]
			#		),
			#		program="RNFtools",
			#		subprogram="bam=>roc",
			#	)

		return vec

	@staticmethod
	def _et_line(readname,vector_of_categories):
		i=0
		intervals=[]
		for j in range(len(vector_of_categories)):
			if vector_of_categories[i]!=vector_of_categories[j]:
				intervals.append("{}:{}-{}".format(vector_of_categories[i],i,j-1))
				i=j
		intervals.append("{}:{}-{}".format(vector_of_categories[i],i,len(vector_of_categories)-1))

		return "{}\t{}".format(readname,",".join(intervals))

	@staticmethod
	def es2et(
				es_fo,
				et_fo,
			):
		"""Convert ES to ET.

		Args:
			es_fo (file): File object for the ES file.
			et_fo (file): File object for the ET file.
		"""

		et_fo.write("# Mapping information for read tuples"+os.linesep)
		et_fo.write("#"+os.linesep)
		et_fo.write("# RN:   read name"+os.linesep)
		et_fo.write("# I:    intervals with asigned categories"+os.linesep)
		et_fo.write("#"+os.linesep)
		et_fo.write("# RN	I"+os.linesep)

		last_rname=""
		for line in es_fo:
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
						voc = Bam._vector_of_categories(single_reads_statistics,rname,nb_of_segments)
						et_fo.write(Bam._et_line(readname=rname,vector_of_categories=voc))
						et_fo.write(os.linesep)

					# nulling
					single_reads_statistics= [
								{
									"U":0,
									"u":0,
									"M":[],
									"m":0,
									"w":0,
									"T":0,
									"t":0,
								} for i in range(rnftools.lavender.MAXIMAL_MAPPING_QUALITY+1)
							]
					last_rname=rname

				####################
				# Unmapped segment #
				####################

				#####
				# U #
				#####
				if category=="U":
					for q in range(len(single_reads_statistics)):
						single_reads_statistics[q]["U"]+=1

				#####
				# u #
				#####
				elif category=="u":
					for q in range(len(single_reads_statistics)):
						single_reads_statistics[q]["u"]+=1

				##################
				# Mapped segment #
				##################

				else:
					mapping_quality=int(mapped.replace("mapped_",""))
					assert 0<=mapping_quality and mapping_quality<=rnftools.lavender.MAXIMAL_MAPPING_QUALITY, mapping_quality

					#####
					# m #
					#####
					if category=="m":
						for q in range(mapping_quality+1):
							single_reads_statistics[q]["m"]+=1
						for q in range(mapping_quality+1,rnftools.lavender.MAXIMAL_MAPPING_QUALITY+1):
							single_reads_statistics[q]["T"]+=1

					#####
					# w #
					#####
					elif category=="w":
						for q in range(mapping_quality+1):
							single_reads_statistics[q]["w"]+=1
						for q in range(mapping_quality+1,rnftools.lavender.MAXIMAL_MAPPING_QUALITY+1):
							single_reads_statistics[q]["t"]+=1

					#####
					# M #
					#####
					else:
						assert category[0]=="M", category
						segment_id=int(category.replace("M_",""))
						for q in range(mapping_quality+1):
							single_reads_statistics[q]["M"].append(segment_id)
						for q in range(mapping_quality+1,rnftools.lavender.MAXIMAL_MAPPING_QUALITY+1):
							single_reads_statistics[q]["t"]+=1

		# last read
		voc = Bam._vector_of_categories(single_reads_statistics,rname,nb_of_segments)
		et_fo.write(Bam._et_line(readname=rname,vector_of_categories=voc))
		et_fo.write(os.linesep)


	def create_et(self):
		"""Create a et file for this BAM file (mapping information about read tuples).

		:raises: ValueError

		"""

		with (gzip.open(self._es_fn,"tr") if self.compress_intermediate_files else open(self._es_fn,"r")) as es_fo:
			with (gzip.open(self._et_fn,"tw+") if self.compress_intermediate_files else open(self._et_fn,"w+")) as et_fo:

				self.es2et(
						es_fo=es_fo,
						et_fo=et_fo,
					)

	############################
	############################
	##
	## ROC
	##
	############################
	############################

	@staticmethod
	def et2roc(
				et_fo,
				roc_fo,
			):
		"""ET to ROC conversion.

		Args:
			et_fo (file): File object for the ET file.
			roc_fo (file): File object for the ROC file.

		:raises: ValueError

		"""

		stats_dicts = [
				{
					"q":q,
					"M":0,
					"w":0,
					"m":0,
					"P":0,
					"U":0,
					"u":0,
					"T":0,
					"t":0,
					"x":0
				}
				for q in range(rnftools.lavender.MAXIMAL_MAPPING_QUALITY+1)
			]

		for line in et_fo:
			line=line.strip()
			if line!="" and line[0]!="#":
				(read_tuple_name,tab,info_categories)=line.partition("\t")
				intervals=info_categories.split(",")
				for interval in intervals:
					category=interval[0]
					(left,colon,right)=interval[2:].partition("-")
					for q in range(int(left),int(right)+1):
						stats_dicts[q][category]+=1

		roc_fo.write("# Numbers of reads in several categories in dependence"+os.linesep)
		roc_fo.write("# on the applied threshold on mapping quality q"+os.linesep)
		roc_fo.write("# "+os.linesep)
		roc_fo.write("# Categories:"+os.linesep)
		roc_fo.write("#        M: Mapped correctly."+os.linesep)
		roc_fo.write("#        w: Mapped to a wrong position."+os.linesep)
		roc_fo.write("#        m: Mapped but should be unmapped."+os.linesep)
		roc_fo.write("#        P: Multimapped."+os.linesep)
		roc_fo.write("#        U: Unmapped and should be unmapped."+os.linesep)
		roc_fo.write("#        u: Unmapped but should be mapped."+os.linesep)
		roc_fo.write("#        T: Thresholded correctly."+os.linesep)
		roc_fo.write("#        t: Thresholded incorrectly."+os.linesep)
		roc_fo.write("#        x: Unknown."+os.linesep)
		roc_fo.write("#"+os.linesep)
		roc_fo.write("# q\tM\tw\tm\tP\tU\tu\tT\tt\tx\tall"+os.linesep)

		l_numbers = []
		for line in stats_dicts:
			numbers = [line["M"],line["w"],line["m"],line["P"],line["U"],line["u"],line["T"],line["t"],line["x"]]
			if numbers != l_numbers:
				roc_fo.write("\t".join(
						[str(line["q"])] + list(map(str,numbers)) + [str(sum(numbers))]
					)+os.linesep)
			l_numbers=numbers

	def create_roc(self):
		"""Create a ROC file for this BAM file.

		:raises: ValueError

		"""

		with (gzip.open(self._et_fn,"tr") if self.compress_intermediate_files else open(self._et_fn,"r")) as et_fo:
			with open(self._roc_fn, "w+") as roc_fo:
				self.et2roc(
						et_fo=et_fo,
						roc_fo=roc_fo,
					)

	############################
	############################
	##
	## GRAPHICS
	##
	############################
	############################

	def create_gp(self):
		"""Create a GnuPlot file for this BAM file."""

		categories_order=[
			("{U}",     "#ee82ee", 'Unmapped correctly'),
			("{u}",     "#ff0000", 'Unmapped incorrectly'),
			("{T}",     "#00ff00", 'Thresholded correctly'),
			("{t}",     "#008800", 'Thresholded incorrectly'),
			("{P}",     "#ffff00", 'Multimapped'),
			("{w}+{x}", "#7f7f7f", 'Mapped, should be unmapped'),
			("{m}",     "#000000", 'Mapped to wrong position'),
			("{M}",     "#0000ff", 'Mapped correctly'),
		]

		plot_lines=[
			'"{roc_fn}" using (( ({x}) )):({y}) lt rgb "{color}" with filledcurve x1 title "{title}", \\'.format(
					roc_fn=self._roc_fn,
					x=rnftools.lavender._format_xxx(self.default_x_axis),
					y=rnftools.lavender._format_xxx('({sum})*100/{{all}}'.format(sum="+".join([c[0] for c in categories_order[i:]]))),
					color=categories_order[i][1],
					title=categories_order[i][2],
				)

			for i in range(len(categories_order))
		]

		plot=os.linesep.join((["plot \\"] + plot_lines))

		with open(self._gp_fn,"w+") as gp:
			gp_content="""
				set title "{{/:Bold=16 {title}}}"

				set x2lab "{x_lab}"
				set log x
				set log x2

				set format x "10^{{%L}}"
				set format x2 "10^{{%L}}"
				set xran  [{xran}]
				set x2ran [{xran}]
				set x2tics
				unset xtics


				set ylab "Part of all reads (%)"

				set format y "%g %%"
				set yran [{yran}]
				set y2ran [{yran}]

				set pointsize 1.5

				set grid ytics lc rgb "#777777" lw 1 lt 0 front
				set grid x2tics lc rgb "#777777" lw 1 lt 0 front

				set datafile separator "\t"
				set palette negative

				set termin svg size {svg_size} enhanced
				set out "{svg_fn}"
				set key spacing 0.8 opaque width -5
				{plot}


				set termin pdf enhanced size {pdf_size} enhanced font 'Arial,12'
				set out "{pdf_fn}"
				set key spacing 0.8 opaque width 0
				{plot}

			""".format(
				svg_fn=self._svg_fn,
				pdf_fn=self._pdf_fn,
				xran="{:.10f}:{:.10f}".format(self.report.default_x_run[0],self.report.default_x_run[1]),
				yran="{:.10f}:{:.10f}".format(self.report.default_y_run[0],self.report.default_y_run[1]),
				svg_size="{},{}".format(self.report.default_svg_size_px[0],self.report.default_svg_size_px[1]),
				pdf_size="{:.10f}cm,{:.10f}cm".format(self.report.default_pdf_size_cm[0],self.report.default_pdf_size_cm[1]),
				title=os.path.basename(self._bam_fn)[:-4],
				plot=plot,
				x_lab=self.default_x_label,
			)
			gp.write(gp_content)


	def create_graphics(self):
		"""Create images related to this BAM file using GnuPlot."""

		smbl.utils.shell('"{}" "{}"'.format(smbl.prog.GNUPLOT5,self._gp_fn))


	############################
	############################
	##
	## HTML
	##
	############################
	############################

	def create_html(self):
		"""Create a HTML page for this BAM file."""

		roc_dicts = []
		with open(self._roc_fn,"r") as roc:
			for line in roc:
				line=line.strip()
				if line!="" and line[0]!="#":
					(q,M,w,m,P,U,u,T,t,x,a)=line.split("\t")
					roc_dict = {
						"q":int(q),
						"M":int(M),
						"w":int(w),
						"m":int(m),
						"P":int(P),
						"U":int(U),
						"u":int(u),
						"T":int(T),
						"t":int(t),
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
						<td>        {P}                         </td>
						<td><small> {P_proc:.2f}        </small></td>
						<td>        {unmapped}                  </td>
						<td><small> {unmapped_proc:.2f} </small></td>
						<td>        {U}                         </td>
						<td><small> {U_proc:.2f}        </small></td>
						<td>        {u}                         </td>
						<td><small> {u_proc:.2f}        </small></td>
						<td>        {T}                         </td>
						<td><small> {T_proc:.2f}        </small></td>
						<td>        {t}                         </td>
						<td><small> {t_proc:.2f}        </small></td>
						<td>        {x}                         </td>
						<td><small> {x_proc:.2f}        </small></td>
						<td>        {sum}                       </td>
						<td>        {prec_proc:.3f}              </td>
					</tr>
				""".format(
						quality       = roc_dict["q"],
						mapped        = roc_dict["M"]+roc_dict["w"]+roc_dict["m"]+roc_dict["P"],
						mapped_proc   = 100.0*(roc_dict["M"]+roc_dict["w"]+roc_dict["m"]+roc_dict["P"])/roc_dict["a"],
						M             = roc_dict["M"],
						M_proc        = 100.0*(roc_dict["M"])/roc_dict["a"],
						w             = roc_dict["w"],
						w_proc        = 100.0*(roc_dict["w"])/roc_dict["a"],
						m             = roc_dict["m"],
						m_proc        = 100.0*(roc_dict["m"])/roc_dict["a"],
						P             = roc_dict["P"],
						P_proc        = 100.0*(roc_dict["P"])/roc_dict["a"],
						unmapped      = roc_dict["U"]+roc_dict["u"]+roc_dict["T"]+roc_dict["t"]+roc_dict["x"],
						unmapped_proc = 100.0*(roc_dict["U"]+roc_dict["u"]+roc_dict["T"]+roc_dict["t"]+roc_dict["x"])/roc_dict["a"],
						U             = roc_dict["U"],
						U_proc        = 100.0*(roc_dict["U"])/roc_dict["a"],
						u             = roc_dict["u"],
						u_proc        = 100.0*(roc_dict["u"])/roc_dict["a"],
						T             = roc_dict["T"],
						T_proc        = 100.0*(roc_dict["T"])/roc_dict["a"],
						t             = roc_dict["t"],
						t_proc        = 100.0*(roc_dict["t"])/roc_dict["a"],
						x             = roc_dict["x"],
						x_proc        = 100.0*(roc_dict["x"])/roc_dict["a"],
						sum           = roc_dict["a"],
						prec_proc      = 100.0*(roc_dict["M"])/(roc_dict["M"]+roc_dict["w"]+roc_dict["m"]+roc_dict["P"]) if (roc_dict["M"]+roc_dict["w"]+roc_dict["m"]+roc_dict["P"]) != 0 else 0,
					)
				for roc_dict in roc_dicts
			])


		with open(self._html_fn,"w+") as html:
			program_info=["No information available (PG header is essing)."]
			for x in smbl.utils.shell(
						'"{samtools}" view -H "{bam}"'.format(
							samtools=smbl.prog.SAMTOOLS,
							bam=self._bam_fn,
						),
						iterable=True
					):
				x=x.strip()

				if x[:3]=="@PG":
					pg_header=x[4:].strip()
					parts=pg_header.split("\t")

					program_info=["<table>"]
					for part in parts:
						(lvalue,colon,rvalue)=part.partition(":")
						program_info.append('<tr><td style="font-weight:bold">{}:</td><td style="text-align:left">{}</td></tr>'.format(lvalue,rvalue))

					program_info.append("</table>")

			html_src="""<!DOCTYPE html>
			<html>
			<head>
				<meta charset="UTF-8" />
				<title>{name}</title>
				<style type="text/css">
					table            {{border-collapse:collapse;}}
					td               {{border: solid #aaaaff 1px;padding:4px;text-align:right;}}
					colgroup, thead  {{border: solid black 2px;padding 2px;}}
					.link_to_top     {{font-size:10pt;}}
					.desc            {{color:#aaa;}}
					.formats         {{text-align:left;margin-bottom:20px;}}
				</style>
			</head>
			<body>
				<h1 id="top">
					{name}
					<span class="link_to_top">
						[<a href="{homepage}">Back to main report</a>]
					</span>
				</h1>

				<p>
					<a href="#info">Information about program</a> -
					<a href="#roctable">ROC table</a> -
					<a href="#graphs">Graphs</a>
				</p>

				<h2 id="info">
					Information about program
					{headline_links}
				</h2>


				{program_info}

				<h2 id="roctable">
					ROC table
					{headline_links}
				</h2>

				<p style="font-size:80%">
					<strong>M</strong>:
						mapped correctly
						<span class="desc">
							(all segments of tuple are mapped once and correctly),
						</span>
					<strong>w</strong>:
						mapped to wrong position
						<span class="desc">
							(at least one segment was mapped to a wrong position),
						</span>
					<strong>m</strong>:
						mapped but should be unmapped
						<span class="desc">
							(at least one segment was mapped but the read should not be mapped),
						</span>
					<strong>P</strong>:
						multimapped
						<span class="desc">
							(read should be mapped but it at least one segment was mapped several
							times and all segments were mapped correctly at least once),
						</span>
					<strong>U</strong>:
						unmapped correctly
						<span class="desc">
							(all segments of the read were correctly marked as unmapped),
						</span>
					<strong>u</strong>:
						unmapped but should be mapped
						<span class="desc">
							(at least one segment was mapped but entire read should be unmapped),
						</span>
					<strong>T</strong>:	
						thresholded correctly
						<span class="desc">
							(read shoud not be mapped),
						</span>
					<strong>t</strong>:
						thresholded incorrectly
						<span class="desc">
							(read should be mapped),
						</span>
					<strong>x</strong>:
						unknown
						<span class="desc">
							(read is probably not reported by mapper)
						</span>
				</p>
				<table>
					<colgroup span="1" style="">
					<colgroup span="2" style="background-color:#ddd">
					<colgroup span="8" style="">
					<colgroup span="2" style="background-color:#ddd">
					<colgroup span="10" style="">
					<colgroup span="1" style="background-color:#ddd">
					<colgroup span="1" style="">
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
							<td>P</td>
							<td>%</td>
							<td>unmapped</td>
							<td>%</td>
							<td>U</td>
							<td>%</td>
							<td>u</td>
							<td>%</td>
							<td>T</td>
							<td>%</td>
							<td>t</td>
							<td>%</td>
							<td>x</td>
							<td>%</td>
							<td>sum</td>
							<td>prec. (%)</td>
						</tr>
					</thead>
					<tbody>
						{tbody}
					</tbody>
				</table>

				<h2 id="graphs">
					Graphs
					{headline_links}
				</h2>

				<div class="formats">
					<img src="{svg}" />
					<br />
					<a href="{pdf}">PDF version</a>
					|
					<a href="{svg}">SVG version</a>
					|
					<a href="{roc}" type="text/csv">ROC file</a>
					|
					<a href="{gp}" type="text/plain">GP file</a>
				</div>


			</body>
			</html>			
			""".format(
					name=self.name,
					tbody=tbody,
					svg=os.path.relpath(
						self._svg_fn,
						os.path.dirname(self._html_fn)
					),
					pdf=os.path.relpath(
						self._pdf_fn,
						os.path.dirname(self._html_fn)
					),
					roc=os.path.relpath(
						self._roc_fn,
						os.path.dirname(self._html_fn)
					),
					gp=os.path.relpath(
						self._gp_fn,
						os.path.dirname(self._html_fn)
					),
					program_info=os.linesep.join(program_info),
					homepage=os.path.relpath(
						self.report.html_fn(),
						os.path.dirname(self._html_fn)
					),
					headline_links='''
						<span class="link_to_top">
							[<a href="#top">Top of this page</a>]
							[<a href="{homepage}">Back to main report</a>]
						</span>
					'''.format(
							homepage=os.path.relpath(
								self.report.html_fn(),
								os.path.dirname(self._html_fn)
							),
						)
				)

			html.write(html_src)
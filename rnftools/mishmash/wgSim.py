import rnftools
from .source import Source

import os
import smbl
import snakemake
import re
import itertools

class WgSim(Source):
	"""Class for the WgSim.

	Single-end reads and pair-end reads simulations are supported. For pair-end simulations,
	ends can have different lengths.
	"""

	#TODO:estimate_unknown_values=False,
	def __init__(self,
			fa,
			coverage=0,
			number_of_reads=0,
			read_length_1=100,
			read_length_2=0,
			other_params="",
			distance=500,
			distance_deviation=50.0,
			rng_seed=1,
			haplotype_mode=False,
		):
		"""
		:param fa: File name of the genome from which reads are created (FASTA file).
		:type  fa: str
		:param coverage: Average coverage of the genome (if number_of_reads specified, then it must be equal to zero).
		:type  coverage: float
		:param number_of_reads: Number of reads (if coverage specified, then it must be equal to zero).
		:type  number_of_reads: int
		:param read_length_1: Length of the first end of a read.
		:type  read_length_1: int
		:param read_length_2: Length of the second end of a read (if zero, then single-end reads are created).
		:type  read_length_2: int
		:param other_params: Other parameters which are used on commandline.
		:type  other_params: str
		:param distance: Mean inner distance between ends.
		:type  distance: int
		:param distance_deviation: Deviation of inner distances between ends.
		:type  distance_deviation: int
		:param rng_seed: Seed for simulator's random number generator.
		:type  rng_seed: int
		:param haplotype_mode: Simulate reads in haplotype mode.
		:type  haplotype_mode: bool
		:raises: ValueError
		"""

		
		if read_length_2==0:
			ends = 1
			#raise NotImplementedError("Since information about strand is not present in WgSim read names, it can be currently used only in pair-end mode.")
		else:
			ends = 2
			self.distance=distance
			self.distance_deviation=distance_deviation

		super().__init__(
				fa=fa,
				ends=ends,
				rng_seed=rng_seed,
			)

		self.read_length_1=read_length_1
		self.read_length_2=read_length_2
		self.haplotype_mode=haplotype_mode
		self.other_params=other_params

		if coverage*number_of_reads!=0:
			raise ValueError("coverage or number_of_reads must be equal to zero")

		self.number_of_reads=number_of_reads
		self.coverage=coverage

		self._tmp_fq1_fn = os.path.join(self._dir,"1.fq")
		if self._ends==2:
			self._tmp_fq2_fn = os.path.join(self._dir,"2.fq")
		else:
			self._tmp_fq2_fn = "/dev/null"


	def get_input(self):
		return [
				smbl.prog.WGSIM,
				self._fa_fn,
				self._fai_fn,
			]

	def get_output(self):
		output = [
				self._tmp_fq1_fn,
				self._fq_fn
			]
		if self._ends==2:
			output.append(self._tmp_fq2_fn)
		return output

	def create_fq(self):
		if self.number_of_reads == 0:
			genome_size=os.stat(self._fa_fn).st_size
			self.number_of_reads=int(self.coverage*genome_size/(self.read_length_1+self.read_length_2))


		if self._ends==2:
			paired_params="-d {dist} -s {dist_dev}".format(
					dist=self.distance,
					dist_dev=self.distance_deviation,
				)
		else:
			paired_params=""

		if self.read_length_2==0:
			fake_read_length_2=42
		else:
			fake_read_length_2=self.read_length_2

		snakemake.shell("""
				{wgsim} \
				-1 {rlen1} \
				-2 {rlen2} \
				-S {rng_seed} \
				-N {nb} \
				{haplotype}\
				{paired_params} \
				{other_params} \
				"{fa}" \
				"{fq1}" \
				"{fq2}" \
				> /dev/null
			""".format(
				wgsim=smbl.prog.WGSIM,
				fa=self._fa_fn,
				fq1=self._tmp_fq1_fn,
				fq2=self._tmp_fq2_fn,
				nb=self.number_of_reads,
				rlen1=self.read_length_1,
				rlen2=fake_read_length_2,
				other_params=self.other_params,
				paired_params=paired_params,
				rng_seed=self._rng_seed,
				haplotype="-h" if self.haplotype_mode else "",
			)
		)
		if self._ends==1:
			self.recode_wgsim_reads(
				old_fq1=self._tmp_fq1_fn,
			)
		else:
			self.recode_wgsim_reads(
				old_fq1=self._tmp_fq1_fn,
				old_fq2=self._tmp_fq2_fn,
			)

	def recode_wgsim_reads(self,old_fq1,old_fq2=None,number_of_reads=10**9):
		wgsim_pattern = re.compile('@(.*)_([0-9]+)_([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9a-f]+)/([12])')
		"""
			WGSIM read name format

		1)  contig name (chromsome name)
		2)  start end 1 (one-based)
		3)  end end 2 (one-based)
		4)  number of errors end 1
		5)  number of substitutions end 1
		6)  number of indels end 1
		5)  number of errors end 2
		6)  number of substitutions end 2
		7)  number of indels end 2
		10) id
		11) pair
		"""
		
		max_seq_len=0

		self.load_fai()
		id_str_size=len(format(number_of_reads,'x'))

		rn_formatter = rnftools.rnfformat.RnFormatter(
				id_str_size=id_str_size,
				source_str_size=2,
				chr_str_size=self.chr_str_size,
				pos_str_size=self.pos_str_size
			)

		# parsing FQ file

		read_id=0
		with open(self._tmp_fq1_fn,"r+") as f_inp_1:
			with open(self._tmp_fq2_fn,"r+") as f_inp_2:
				with open(self._fq_fn,"w+") as f_out:

					if self._ends==2:
						f_inp_2 = open(self._tmp_fq2_fn,"r+")
						read_suffix="[pair-end,wgsim]"
					else:
						read_suffix="[single-end,wgsim]"
					i=[0,0]
					(fq_line1,fq_line2,fq_line3,fq_line4)=(["",""],["",""],["",""],["",""])
					for line_a in f_inp_1:
						line=[line_a]
						if self._ends==2:
							line_b=f_inp_2.readline()
							line.append(line_b)

						for end_id in range(self._ends):
							if i[end_id]%4==0:
								m = wgsim_pattern.search(line[end_id])
								if m is None:
									raise ValueError("Read '{}' was not generated by WgSim.".format(line[end_id][1:]))

								contig_name     = m.group(1)
								start_1         = int(m.group(2))
								end_2           = int(m.group(3))
								errors_1        = int(m.group(4))
								substitutions_1 = int(m.group(5))
								indels_1        = int(m.group(6))
								errors_2        = int(m.group(7))
								substitutions_2 = int(m.group(8))
								indels_2        = int(m.group(9))
								read_id_w       = int(m.group(10),16)
								pair            = int(m.group(11))

								chr_id = self.dict_chr_ids[contig_name] if self.dict_chr_ids!={} else "0"

							elif i[end_id]%4==1:
								fq_line2[end_id]=line[end_id]
								read_length=len(fq_line2[end_id].strip())

								if start_1<end_2:
									direction_1="F"
									direction_2="R"
								else:
									direction_1="R"
									direction_2="F"

								segments=[
									rnftools.rnfformat.Segment(
										source=self.source_id,
										chr=chr_id,
										direction=direction_1,
										left=start_1,
										right=0),

									rnftools.rnfformat.Segment(
										source=self.source_id,
										chr=chr_id,
										direction=direction_2,
										left=0,
										right=end_2)
								]

								read = rnftools.rnfformat.Read(segments=segments,read_id=read_id+1,suffix=read_suffix)
								new_read_name = rn_formatter.process_read(read)

								if self._ends==1:
									ends_suffix=""
									read_id+=1 
								else:
									if end_id==0:
										ends_suffix="/1"
									else:
										ends_suffix="/2"
										read_id+=1

								fq_line1[end_id]="".join(["@",new_read_name,ends_suffix,os.linesep])

							elif i[end_id]%4==2:
								fq_line3[end_id]=line[end_id]

							else:
								fq_line4[end_id]=line[end_id]

								f_out.write(fq_line1[end_id])
								f_out.write(fq_line2[end_id])
								f_out.write(fq_line3[end_id])
								f_out.write(fq_line4[end_id])

							i[end_id]+=1					


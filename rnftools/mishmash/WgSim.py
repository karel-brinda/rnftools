import rnftools
from .Source import *

import os
import snakemake
import re
import itertools


class WgSim(Source):
	"""Class for the WGsim (https://github.com/lh3/wgsim).

	Single-end and pair-end simulations are supported. For pair-end simulations,
	reads can have different lengths.

	Args:
		fasta (str): File name of the genome from which reads are created (FASTA file).
		sequences (set of int or str): FASTA sequences to extract. Sequences can be specified either by their ids, or by their names.
		coverage (float): Average coverage of the genome (if number_of_read_tuples specified, then it must be equal to zero).
		number_of_read_tuples (int): Number of read tuples (if coverage specified, then it must be equal to zero). Corresponding WGsim parameter: ``-N``.
		read_length_1 (int): Length of the first read. Corresponding WGsim parameter: ``-1``.
		read_length_2 (int): Length of the second read (if zero, then single-end reads are simulated). Corresponding WGsim parameter: ``-2``.
		distance (int): Mean outer distance of reads. Corresponding WGsim parameter: ``-d``.
		distance_deviation (int): Standard deviation of outer distances of reads. Corresponding WGsim parameter: ``-s``.
		rng_seed (int): Seed for simulator's random number generator. Corresponding WGsim parameter: ``-S``.
		haploid_mode (bools): Simulate reads in haploid mode. Corresponding WGsim parameter: ``-h``.
		error_rate (float): Sequencing error rate (sequencing errors). Corresponding WGsim parameter: ``-e``.
		mutation_rate (float): Mutation rate. Corresponding WGsim parameter: ``-r``.
		indels (float): Rate of indels in mutations. Corresponding WGsim parameter: ``-R``.
		prob_indel_ext (float): Probability that an indel is extended. Corresponding WGsim parameter: ``-X``.
		other_params (str): Other parameters on commandline.

	Raises:
		ValueError
	"""

	# TODO:estimate_unknown_values=False,
	def __init__(
			self,
			fasta,
			sequences=None,
			coverage=0,
			number_of_read_tuples=0,
			read_length_1=100,
			read_length_2=0,
			distance=500,
			distance_deviation=50.0,
			rng_seed=1,
			haploid_mode=False,
			error_rate=0.020,
			mutation_rate=0.001,
			indels=0.15,
			prob_indel_ext=0.3,
			other_params="",
	):

		if read_length_2 == 0:
			ends = 1
		else:
			ends = 2
			self.distance = distance
			self.distance_deviation = distance_deviation

		super().__init__(
			fasta=fasta,
			sequences=sequences,
			reads_in_tuple=ends,
			rng_seed=rng_seed,
		)

		self.read_length_1 = read_length_1
		self.read_length_2 = read_length_2
		self.haploid_mode = haploid_mode
		self.other_params = other_params

		self.error_rate = error_rate
		self.mutation_rate = mutation_rate
		self.indels = indels
		self.prob_indel_ext = prob_indel_ext

		if coverage * number_of_read_tuples != 0:
			rnftools.utils.error("coverage or number_of_read_tuples must be equal to zero", program="RNFtools",
				subprogram="MIShmash", exception=ValueError)

		self.number_of_read_tuples = number_of_read_tuples
		self.coverage = coverage

		self._tmp_fq1_fn = os.path.join(self._dir, "wgsim_files.1.fq")
		if self._reads_in_tuple == 2:
			self._tmp_fq2_fn = os.path.join(self._dir, "wgsim_files.2.fq")
		else:
			self._tmp_fq2_fn = "/dev/null"

	def get_input(self):
		return [
			self._fa_fn,
			self._fai_fn,
		]

	def get_output(self):
		output = [
			self._tmp_fq1_fn,
			self._fq_fn
		]
		if self._reads_in_tuple == 2:
			output.append(self._tmp_fq2_fn)
		return output

	def create_fq(self):
		if self.number_of_read_tuples == 0:
			genome_size = os.stat(self._fa_fn).st_size
			self.number_of_read_tuples = int(self.coverage * genome_size / (self.read_length_1 + self.read_length_2))

		if self._reads_in_tuple == 2:
			paired_params = "-d {dist} -s {dist_dev}".format(
				dist=self.distance,
				dist_dev=self.distance_deviation,
			)
		else:
			paired_params = ""

		if self.read_length_2 == 0:
			fake_read_length_2 = 42
		else:
			fake_read_length_2 = self.read_length_2

		rnftools.utils.shell("""
				{wgsim} \
				-1 {rlen1} \
				-2 {rlen2} \
				-S {rng_seed} \
				-N {nb} \
				-e {error_rate} \
				-r {mutation_rate} \
				-R {indels} \
				-X {prob_indel_ext} \
				{haploid}\
				{paired_params} \
				{other_params} \
				"{fa}" \
				"{fq1}" \
				"{fq2}" \
				> /dev/null
			""".format(
			wgsim="wgsim",
			fa=self._fa_fn,
			fq1=self._tmp_fq1_fn,
			fq2=self._tmp_fq2_fn,
			nb=self.number_of_read_tuples,
			rlen1=self.read_length_1,
			rlen2=fake_read_length_2,
			other_params=self.other_params,
			paired_params=paired_params,
			rng_seed=self._rng_seed,
			haploid="-h" if self.haploid_mode else "",
			error_rate=self.error_rate,
			mutation_rate=self.mutation_rate,
			indels=self.indels,
			prob_indel_ext=self.prob_indel_ext,
		)
		)

		with open(self._fai_fn) as fai_fo:
			with open(self._fq_fn, "w+") as fq_fo:
				self.recode_wgsim_reads(
					rnf_fastq_fo=fq_fo,
					fai_fo=fai_fo,
					genome_id=self.genome_id,
					wgsim_fastq_1_fn=self._tmp_fq1_fn,
					wgsim_fastq_2_fn=self._tmp_fq2_fn if self._reads_in_tuple == 2 else None,
					number_of_read_tuples=10 ** 9,
				)

	@staticmethod
	def recode_wgsim_reads(
			rnf_fastq_fo,
			fai_fo,
			genome_id,
			wgsim_fastq_1_fn,
			wgsim_fastq_2_fn=None,
			number_of_read_tuples=10 ** 9,
	):
		"""Convert WgSim FASTQ files to RNF FASTQ files.

		Args:
			rnf_fastq_fo (file): File object of the target RNF file.
			fai_fo (file): File object of FAI index of the reference genome.
			genome_id (int): RNF genome ID.
			wgsim_fastq_1_fn (str): File name of the first WgSim FASTQ file.
			wgsim_fastq_2_fn (str): File name of the second WgSim FASTQ file.
			number_of_read_tuples (int): Expected number of read tuples (to estimate widths).
		"""

		wgsim_pattern = re.compile(
			'@(.*)_([0-9]+)_([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9a-f]+)/([12])')
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

		fai_index = rnftools.utils.FaIdx(fai_fo)
		read_tuple_id_width = len(format(number_of_read_tuples, 'x'))

		last_read_tuple_name = None

		fq_creator = rnftools.rnfformat.FqCreator(
			fastq_fo=rnf_fastq_fo,
			read_tuple_id_width=read_tuple_id_width,
			genome_id_width=2,
			chr_id_width=fai_index.chr_id_width,
			coor_width=fai_index.coor_width,
			info_reads_in_tuple=True,
			info_simulator="wgsim",
		)

		reads_in_tuple = 2
		if wgsim_fastq_2_fn is None:
			reads_in_tuple = 1

		i = 0
		with open(wgsim_fastq_1_fn, "r+") as f_inp_1:
			if reads_in_tuple == 2:
				# todo: close file
				f_inp_2 = open(wgsim_fastq_2_fn)

			for line_a in f_inp_1:
				lines = [line_a.strip()]
				if reads_in_tuple == 2:
					lines.append(f_inp_2.readline().strip())

				if i % 4 == 0:
					segments = []
					# bases=[]
					# qualities=[]

					m = wgsim_pattern.search(lines[0])
					if m is None:
						rnftools.utils.error("Read tuple '{}' was not generated by WgSim.".format(
							lines[0][1:]),
							program="RNFtools",
							subprogram="MIShmash",
							exception=ValueError
						)

					contig_name = m.group(1)
					start_1 = int(m.group(2))
					end_2 = int(m.group(3))
					errors_1 = int(m.group(4))
					substitutions_1 = int(m.group(5))
					indels_1 = int(m.group(6))
					errors_2 = int(m.group(7))
					substitutions_2 = int(m.group(8))
					indels_2 = int(m.group(9))
					read_tuple_id_w = int(m.group(10), 16)
					pair = int(m.group(11))

					chr_id = fai_index.dict_chr_ids[contig_name] if fai_index.dict_chr_ids != {} else "0"

					if start_1 < end_2:
						direction_1 = "F"
						direction_2 = "R"
					else:
						direction_1 = "R"
						direction_2 = "F"

					segment1 = rnftools.rnfformat.Segment(
						genome_id=genome_id,
						chr_id=chr_id,
						direction=direction_1,
						left=start_1,
						right=0,
					)

					segment2 = rnftools.rnfformat.Segment(
						genome_id=genome_id,
						chr_id=chr_id,
						direction=direction_2,
						left=0,
						right=end_2,
					)

				elif i % 4 == 1:
					bases = lines[0]
					if reads_in_tuple == 2:
						bases2 = lines[1]

				elif i % 4 == 2:
					pass

				elif i % 4 == 3:
					qualities = lines[0]
					if reads_in_tuple == 2:
						qualities2 = lines[1]

					if reads_in_tuple == 1:
						fq_creator.add_read(
							read_tuple_id=i // 4 + 1,
							bases=bases,
							qualities=qualities,
							segments=[segment1, segment2],
						)
					else:
						fq_creator.add_read(
							read_tuple_id=i // 4 + 1,
							bases=bases,
							qualities=qualities,
							segments=[segment1],
						)
						fq_creator.add_read(
							read_tuple_id=i // 4 + 1,
							bases=bases2,
							qualities=qualities2,
							segments=[segment2],
						)

				i += 1

		fq_creator.flush_read_tuple()

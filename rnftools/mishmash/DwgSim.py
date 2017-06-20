import rnftools
from .Source import *

import os
import snakemake
import re


class DwgSim(Source):
	"""Class for DWGsim (https://github.com/nh13/DWGSIM/wiki).

	Both single-end and paired-end simulations are supported. In paired-end simulations,
	reads can have different lengths. Note that there is a bug in DWGsim documentation:
	coordinates are 1-based.

	Args:
		fasta (str): File name of the genome from which reads are created (FASTA file). 
		sequences (set of int or str): FASTA sequences to extract. Sequences can be specified either by their ids, or by their names.
		coverage (float): Average coverage of the genome (if number_of_reads specified,
			then it must be equal to zero).
			Corresponding DWGsim parameter: ``-C``.
		number_of_read_tuples (int): Number of read tuples (if coverage specified, then
			it must be equal to zero).
			Corresponding DWGsim parameter: ``-N``.
		read_length_1 (int): Length of the first read.
			Corresponding DWGsim parameter: ``-1``.
		read_length_2 (int): Length of the second read (if zero, then single-end
			simulation performed).
			Corresponding DWGsim parameter: ``-2``.
		distance (int): Mean inner distance between reads.
			Corresponding DWGsim parameter: ``-d``.
		distance_deviation (int): Standard deviation of inner distances between both reads.
			Corresponding DWGsim parameter: ``-s``.
		rng_seed (int): Seed for simulator's random number generator.
			Corresponding DWGsim parameter: ``-z``.
		haploid_mode (bools): Simulate reads in haploid mode.
			Corresponding DWGsim parameter: ``-H``.
		error_rate_1 (float): Sequencing error rate in the first read.
			Corresponding DWGsim parameter: ``-e``.
		error_rate_2 (float): Sequencing error rate in the second read.
			Corresponding DWGsim parameter: ``-E``.
		mutation_rate (float): Mutation rate.
			Corresponding DWGsim parameter: ``-e``.
		indels (float): Rate of indels in mutations.
			Corresponding DWGsim parameter: ``-R``.
		prob_indel_ext (float): Probability that an indel is extended.
			Corresponding DWGsim parameter: ``-X``.
		estimate_unknown_values (bool): Estimate unknown values (coordinates missing in
			DWGsim output).
		other_params (str): Other parameters which are used on command-line.

	Raises:
		ValueError
	"""

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
			error_rate_1=0.020,
			error_rate_2=0.020,
			mutation_rate=0.001,
			indels=0.15,
			prob_indel_ext=0.3,
			estimate_unknown_values=False,
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
		self.other_params = other_params

		if coverage * number_of_read_tuples != 0:
			rnftools.utils.error("coverage or number_of_read_tuples must be equal to zero", program="RNFtools",
				subprogram="MIShmash", exception=ValueError)

		self.number_of_read_tuples = number_of_read_tuples
		self.coverage = coverage

		self.haploid_mode = haploid_mode
		self.error_rate_1 = error_rate_1
		self.error_rate_2 = error_rate_2
		self.mutation_rate = mutation_rate
		self.indels = indels
		self.prob_indel_ext = prob_indel_ext

		self.estimate_unknown_values = estimate_unknown_values

		self.dwg_prefix = os.path.join(
			self.get_dir(),
			"dwgsim_files.{}.{}".format("se" if self.number_of_read_tuples == 1 else "pe", self.genome_id)
		)

	def get_input(self):
		return [
			self._fa_fn,
			self._fai_fn,
		]

	def get_output(self):
		return [
			self.dwg_prefix + ".bwa.read1.fastq",
			self.dwg_prefix + ".bwa.read2.fastq",
			self.dwg_prefix + ".bfast.fastq",
			self.dwg_prefix + ".mutations.vcf",
			self.dwg_prefix + ".mutations.txt",
			self._fq_fn,
		]

	def create_fq(self):
		if self.coverage == 0 and self.number_of_read_tuples == 0:
			for x in self.get_output():
				with open(x, "w+") as f:
					f.write(os.linesep)

		else:

			if self.number_of_read_tuples == 0:
				genome_size = os.stat(self._fa_fn).st_size
				self.number_of_read_tuples = int(
					self.coverage * genome_size / (self.read_length_1 + self.read_length_2))

			if self._reads_in_tuple == 2:
				paired_params = "-d {dist} -s {dist_dev}".format(
					dist=self.distance,
					dist_dev=self.distance_deviation,
				)
			else:
				paired_params = ""

			rnftools.utils.shell("""
					"{dwgsim}" \
					-1 {rlen1} \
					-2 {rlen2} \
					-z {rng_seed} \
					-y 0 \
					-N {nb} \
					-e {error_rate_1} \
					-E {error_rate_2} \
					-r {mutation_rate} \
					-R {indels} \
					-X {prob_indel_ext} \
					{haploid} \
					{paired_params} \
					{other_params} \
					"{fa}" \
					"{pref}" \
					> /dev/null
				""".format(
				dwgsim="dwgsim",
				fa=self._fa_fn,
				pref=self.dwg_prefix,
				nb=self.number_of_read_tuples,
				rlen1=self.read_length_1,
				rlen2=self.read_length_2,
				other_params=self.other_params,
				paired_params=paired_params,
				rng_seed=self._rng_seed,
				haploid="-H" if self.haploid_mode else "",
				error_rate_1=self.error_rate_1,
				error_rate_2=self.error_rate_2,
				mutation_rate=self.mutation_rate,
				indels=self.indels,
				prob_indel_ext=self.prob_indel_ext,
			)
			)
			with open(self._fq_fn, "w+") as fastq_fo:
				with open(self._fai_fn) as fai_fo:
					self.recode_dwgsim_reads(
						dwgsim_prefix=self.dwg_prefix,
						fastq_rnf_fo=fastq_fo,
						fai_fo=fai_fo,
						genome_id=self.genome_id,
						number_of_read_tuples=10 ** 9,
						# allow_unmapped=False,
						estimate_unknown_values=self.estimate_unknown_values,
					)

	@staticmethod
	def recode_dwgsim_reads(
			dwgsim_prefix,
			fastq_rnf_fo,
			fai_fo,
			genome_id,
			estimate_unknown_values,
			number_of_read_tuples=10 ** 9,
	):
		"""Convert DwgSim FASTQ file to RNF FASTQ file.

		Args:
			dwgsim_prefix (str): DwgSim prefix of the simulation (see its commandline parameters).
			fastq_rnf_fo (file): File object of RNF FASTQ.
			fai_fo (file): File object for FAI file of the reference genome.
			genome_id (int): RNF genome ID to be used.
			estimate_unknown_values (bool): Estimate unknown values (right coordinate of each end).
			number_of_read_tuples (int): Estimate of number of simulated read tuples (to set width).
		"""

		dwgsim_pattern = re.compile(
			'@(.*)_([0-9]+)_([0-9]+)_([01])_([01])_([01])_([01])_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_(([0-9abcdef])+)')

		###
		# DWGSIM read name format
		# 		
		# 1)  contig name (chromsome name)
		# 2)  start end 1 (one-based)
		# 3)  start end 2 (one-based)
		# 4)  strand end 1 (0 - forward, 1 - reverse)
		# 5)  strand end 2 (0 - forward, 1 - reverse)
		# 6)  random read end 1 (0 - from the mutated reference, 1 - random)
		# 7)  random read end 2 (0 - from the mutated reference, 1 - random)
		# 8)  number of sequencing errors end 1 (color errors for colorspace)
		# 9)  number of SNPs end 1
		# 10) number of indels end 1
		# 11) number of sequencing errors end 2 (color errors for colorspace)
		# 12) number of SNPs end 2
		# 13) number of indels end 2
		# 14) read number (unique within a given contig/chromosome)
		###

		fai_index = rnftools.utils.FaIdx(fai_fo=fai_fo)
		read_tuple_id_width = len(format(number_of_read_tuples, 'x'))

		# parsing FQ file
		read_tuple_id = 0
		last_read_tuple_name = None
		old_fq = "{}.bfast.fastq".format(dwgsim_prefix)

		fq_creator = rnftools.rnfformat.FqCreator(
			fastq_fo=fastq_rnf_fo,
			read_tuple_id_width=read_tuple_id_width,
			genome_id_width=2,
			chr_id_width=fai_index.chr_id_width,
			coor_width=fai_index.coor_width,
			info_reads_in_tuple=True,
			info_simulator="dwgsim",
		)

		i = 0
		with open(old_fq, "r+") as f1:
			for line in f1:
				if i % 4 == 0:
					read_tuple_name = line[1:].strip()
					if read_tuple_name != last_read_tuple_name:
						new_tuple = True
						if last_read_tuple_name is not None:
							read_tuple_id += 1
					else:
						new_tuple = False

					last_read_tuple_name = read_tuple_name
					m = dwgsim_pattern.search(line)
					if m is None:
						rnftools.utils.error(
							"Read tuple '{}' was not created by DwgSim.".format(line[1:]),
							program="RNFtools",
							subprogram="MIShmash",
							exception=ValueError,
						)

					contig_name = m.group(1)
					start_1 = int(m.group(2))
					start_2 = int(m.group(3))
					direction_1 = "F" if int(m.group(4)) == 0 else "R"
					direction_2 = "F" if int(m.group(5)) == 0 else "R"
					# random_1 = bool(m.group(6))
					# random_2 = bool(m.group(7))
					# seq_err_1 = int(m.group(8))
					# snp_1 = int(m.group(9))
					# indels_1 = int(m.group(10))
					# seq_err_2 = int(m.group(11))
					# snp_2 = int(m.group(12))
					# indels_2 = int(m.group(13))
					# read_tuple_id_dwg = int(m.group(14), 16)

					chr_id = fai_index.dict_chr_ids[contig_name] if fai_index.dict_chr_ids != {} else "0"

				elif i % 4 == 1:
					bases = line.strip()

					if new_tuple:

						segment = rnftools.rnfformat.Segment(
							genome_id=genome_id,
							chr_id=chr_id,
							direction=direction_1,
							left=start_1,
							right=start_1 + len(bases) - 1 if estimate_unknown_values else 0,
						)

					else:

						segment = rnftools.rnfformat.Segment(
							genome_id=genome_id,
							chr_id=chr_id,
							direction=direction_2,
							left=start_2,
							right=start_2 + len(bases) - 1 if estimate_unknown_values else 0,
						)

				elif i % 4 == 2:
					pass

				elif i % 4 == 3:
					qualities = line.strip()
					fq_creator.add_read(
						read_tuple_id=read_tuple_id,
						bases=bases,
						qualities=qualities,
						segments=[segment],
					)

				i += 1

		fq_creator.flush_read_tuple()

import rnftools
from .Source import Source

import snakemake
import os


class MasonIllumina(Source):
	"""Class for the Mason  - Illumina mode (https://www.seqan.de/projects/mason/).

	Single-end reads and pair-end reads simulations are supported. For pair-end simulations,
	lengths of both ends must be equal.

	Args:
		fasta (str): File name of the genome from which read tuples are created (FASTA file). Corresponding Mason parameter: ``-ir, --input-reference``.
		sequences (set of int or str): FASTA sequences to extract. Sequences can be specified either by their ids, or by their names.
		coverage (float): Average coverage of the genome (if number_of_reads specified, then it must be equal to zero).
		number_of_read_tuples (int): Number of read tuples (if coverage specified, then it must be equal to zero). Corresponding Mason parameter: ``-n, --num-fragments``.
		read_length_1 (int): Length of the first read. Corresponding Mason parameter: ``--illumina-read-length``.
		read_length_2 (int): Length of the read read (if zero, then single-end reads are simulated). Corresponding Mason parameter: ``--illumina-read-length``.
		distance (int): Mean inner distance between reads. Corresponding Mason parameter: ``--fragment-mean-size``.
		distance_deviation (int): Standard devation of inner distances between reads. Corresponding Mason parameter: ``--fragment-size-std-dev``.
		rng_seed (int): Seed for simulator's random number generator. Corresponding Mason parameter: ``--seed``.
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
			distance_deviation=50,
			rng_seed=1,
			other_params="",
	):

		if read_length_2 == 0:
			ends = 1
		else:
			ends = 2
			self.distance = int(distance)
			self.distance_deviation = int(distance_deviation)
			if read_length_1 != read_length_2:
				rnftools.utils.error(
					"Mason can simulate only pairs with equal lengths",
					program="RNFtools",
					subprogram="MIShmash",
					exception=ValueError,
				)

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
			rnftools.utils.error(
				"coverage or number_of_read_tuples must be equal to zero",
				program="RNFtools",
				subprogram="MIShmash",
				exception=ValueError
			)

		self.number_of_read_tuples = number_of_read_tuples
		self.coverage = coverage

		self.mason_prefix = os.path.join(
			self.get_dir(),
			"mason_files.{}.{}".format("se" if self.number_of_read_tuples == 1 else "pe", self.genome_id)
		)

		self._sam_fn = self.mason_prefix + ".sam"

	def get_input(self):
		return [
			self._fa_fn,
			self._fai_fn,
		]

	def get_output(self):
		if self._reads_in_tuple == 1:
			fqs = [self.mason_prefix + "1.fq"]
		else:
			fqs = [self.mason_prefix + "1.fq", self.mason_prefix + "2.fq"]
		return [self._fq_fn,
				   self._sam_fn,
			   ] + fqs

	def create_fq(self):
		if self.coverage == 0 and self.number_of_read_tuples == 0:
			for x in self.get_output():
				with open(x, "w+") as f:
					f.write(os.linesep)

		else:
			if self.coverage == 0:
				genome_size = os.stat(self._fa_fn).st_size
				self.coverage = 1.0 * self.number_of_read_tuples * (self.read_length_1 + self.read_length_2) / (
					0.8 * genome_size)

			if self._reads_in_tuple == 2:
				paired_params = '--fragment-mean-size {dist} --fragment-size-std-dev {dist_dev} -or "{fq2}"'.format(
					dist=self.distance,
					dist_dev=self.distance_deviation,
					fq2=self.mason_prefix + "2.fq",
				)
			else:
				paired_params = ""

			command = """
					"{mason}" \
						-n {number_of_read_tuples} \
						-ir "{fasta}" \
						--illumina-read-length {rlen} \
						--seed {rng_seed} \
						-o "{fq1}" \
						-oa "{sam}" \
						{paired_params} \
						{other_params} \
						> /dev/null
				""".format(
				mason="mason_simulator",
				paired_params=paired_params,
				fasta=self._fa_fn,
				rlen=self.read_length_1,
				other_params=self.other_params,
				number_of_read_tuples=self.number_of_read_tuples,
				fq1=self.mason_prefix + "1.fq",
				rng_seed=self._rng_seed,
				sam=self._sam_fn,
			)

			rnftools.utils.shell(command)

			with open(self._fq_fn, "w+") as fq_fo:
				with open(self._fai_fn) as fai_fo:
					self.recode_sam_reads(
						sam_fn=self._sam_fn,
						fastq_rnf_fo=fq_fo,
						fai_fo=fai_fo,
						genome_id=self.genome_id,
						number_of_read_tuples=10 ** 9,
						simulator_name="mason",
						allow_unmapped=False,
					)

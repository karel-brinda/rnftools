import rnftools
from .Source import Source

import smbl
import snakemake
import os

class MasonIllumina(Source):
	"""Class for the Mason (Illumina mode).

	Single-end reads and pair-end reads simulations are supported. For pair-end simulations,
	lengths of both ends must be equal.

	Args:
		fasta (str): File name of the genome from which read tuples are created (FASTA file).
		coverage (float): Average coverage of the genome.
		number_of_read_tuples (int): Number of read tuples.
		read_length_1 (int): Length of the first end of a read tuple.
		read_length_2 (int): Length of the second end of a read tuple (if zero, then single-end reads are created).
		other_params (str): Other parameters which are used on commandline.
		distance (int): Mean inner distance between ends.
		distance_deviation (int): Devation of inner distances between ends.
		rng_seed (int): Seed for simulator's random number generator.

	Raises:
		ValueError
	"""

	def __init__(self,
				fasta,
				coverage=0,
				number_of_read_tuples=0,
				read_length_1=100,
				read_length_2=0,
				other_params="",
				distance=500,
				distance_deviation=50,
				rng_seed=1
			):

		if read_length_2==0:
			ends = 1
		else:
			ends = 2
			self.distance=int(distance)
			self.distance_deviation=int(distance_deviation)
			if read_length_1!=read_length_2:
				smbl.messages.error("mason can simulate only pairs with equal lengths",program="RNFtools",subprogram="MIShmash",exception=ValueError)
		
		super().__init__(
				fasta=fasta,
				reads_in_tuple=ends,
				rng_seed=rng_seed,
			)

		self.read_length_1=read_length_1
		self.read_length_2=read_length_2
		self.other_params=other_params


		if coverage*number_of_read_tuples!=0:
			smbl.messages.error("coverage or number_of_read_tuples must be equal to zero",program="RNFtools",subprogram="MIShmash",exception=ValueError)

		self.number_of_read_tuples=number_of_read_tuples
		self.coverage=coverage

		self.mason_prefix=os.path.join(
			self.get_dir(),
			"tmp.{}".format(self.genome_id)
		)

		self._sam_fn = self.mason_prefix+".sam"

	def get_input(self):
		return [
				smbl.prog.MASON_SIMULATOR,
				smbl.prog.SAMTOOLS,
				self._fa_fn,
				self._fai_fn,
			]

	def get_output(self):
		return 	[
				self._fq_fn,
				self._sam_fn,
				self.mason_prefix+"1.fq" if self._reads_in_tuple==1 else
					[self.mason_prefix+"1.fq",self.mason_prefix+"2.fq"],
			]


	def create_fq(self):
		if self.coverage == 0:
			genome_size=os.stat(self._fa_fn).st_size
			self.coverage = 1.0 * self.number_of_read_tuples * (self.read_length_1+self.read_length_2) / (0.8 * genome_size)

		if self._reads_in_tuple==2:
			paired_params='--fragment-mean-size {dist} --fragment-size-std-dev {dist_dev} -or "{fq2}"'.format(
					dist=self.distance,
					dist_dev=self.distance_deviation,
					fq2=self.mason_prefix+"2.fq",
				)
		else:
			paired_params=""

		command ="""
				{mason} \
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
				mason=smbl.prog.MASON_SIMULATOR,
				paired_params=paired_params,
				fasta=self._fa_fn,
				rlen=self.read_length_1,
				other_params=self.other_params,
				number_of_read_tuples=self.number_of_read_tuples,
				fq1=self.mason_prefix+"1.fq",
				rng_seed=self._rng_seed,
				sam=self._sam_fn,
			)

		smbl.utils.shell(command)

		with open(self._fq_fn,"w+") as fq_fo:
			with open(self._fai_fn) as fai_fo:
				self.recode_sam_reads(
					sam_fn=self._sam_fn,
					fastq_rnf_fo=fq_fo,
					fai_fo=fai_fo,
					genome_id=self.genome_id,
					number_of_read_tuples=10**9,
					simulator_name="mason",
					allow_unmapped=False,
				)

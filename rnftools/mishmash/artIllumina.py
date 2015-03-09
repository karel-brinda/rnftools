import rnftools
from .source import Source

import smbl
import snakemake
import os

#
# AUXILIARY FUNCTIONS
#

class ArtIllumina(Source):
	"""Class for the ART Illumina.

	Single-end reads and pair-end reads simulations are supported. For pair-end simulations,
	lengths of both ends must be equal.
	"""

	def __init__(self,
			fa,
			coverage=0,
			number_of_reads=0,
			read_length_1=100,
			read_length_2=0,
			other_params="",
			distance=500,
			distance_deviation=50.0,
			rng_seed=1
		):
		"""
		:param fa: File name of the genome from which reads are created (FASTA file).
		:type  fa: str
		:param coverage: Average coverage of the genome.
		:type  coverage: float
		:param read_length_1: Length of the first end of a read.
		:type  read_length_1: int
		:param read_length_2: Length of the second end of a read (if zero, then single-end reads are created).
		:type  read_length_2: int
		:param other_params: Other parameters which are used on commandline.
		:type  other_params: str
		:param distance: Mean inner distance between ends.
		:type  distance: int
		:param distance_deviation: Devation of inner distances between ends.
		:type  distance_deviation: int
		:param rng_seed: Seed for simulator's random number generator.
		:type  rng_seed: int
		:raises: ValueError
		"""

		if read_length_2==0:
			ends = 1
		else:
			ends = 2
			self.distance=distance
			self.distance_deviation=distance_deviation
			if read_length_1!=read_length_2:
				raise ValueError("art_illumina can simulate only pairs with equal lengths")

		
		super().__init__(
				fa=fa,
				ends=ends,
				rng_seed=rng_seed,
			)

		self.read_length_1=read_length_1
		self.read_length_2=read_length_2
		self.other_params=other_params


		if coverage*number_of_reads!=0:
			raise ValueError("coverage or number_of_reads must be equal to zero")

		self.number_of_reads=number_of_reads
		self.coverage=coverage

		self.art_prefix=os.path.join(
			self.get_dir(),
			"tmp.{}".format(self.source_id)
		)

		self._sam1_fn = self.art_prefix+".sam"
		self._sam2_fn = self.art_prefix+".corrected.sam"

	def get_input(self):
		return [
				smbl.prog.ART_ILLUMINA,
				smbl.prog.SAMTOOLS,
				self._fa_fn,
				self._fai_fn,
			]

	def get_output(self):
		return 	[
				self._fq_fn,
				self._sam1_fn,
				self._sam2_fn,
				self.art_prefix+".fq" if self._ends==1 else
					[self.art_prefix+"1.fq",self.art_prefix+"2.fq"],
			]


	def create_fq(self):
		"""Perform the simulation."""
		if self.coverage == 0:
			genome_size=os.stat(self._fa_fn).st_size
			self.coverage = 1.0 * self.number_of_reads * (self.read_length_1+self.read_length_2) / (0.8 * genome_size)

		if self._ends==2:
			paired_params="-p -m {dist} -s {dist_dev}".format(
					dist=self.distance,
					dist_dev=self.distance_deviation,
				)
		else:
			paired_params=""

		command_1 ="""
				{art_il} -sam -na \
					-i "{fasta}" \
					-l {rlen} \
					-rs {rng_seed} \
					-f {coverage} \
					-o "{o_pref}" \
					{paired_params} \
					{other_params} \
					> /dev/null
			""".format(
				art_il=smbl.prog.ART_ILLUMINA,
				paired_params=paired_params,
				fasta=self._fa_fn,
				rlen=self.read_length_1,
				other_params=self.other_params,
				coverage=self.coverage,
				o_pref=self.art_prefix,
				rng_seed=self._rng_seed,
			)

		# correction of header (bug in ART)
		command_2 ="""
			cat "{sam_1}" | \
			grep -v ^@ | \
			"{samtools}" view -h -T "{fa}" - \
			> "{sam_2}"
		""".format(
				samtools=smbl.prog.SAMTOOLS,
				sam_1=self._sam1_fn,
				sam_2=self._sam2_fn,
				fa=self._fa_fn,
		)

		snakemake.shell(command_1)
		snakemake.shell(command_2)

		self.recode_sam_reads(
			sam=self._sam2_fn,
			simulator_name="art-illumina",
		)

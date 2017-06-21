import rnftools
from .Source import Source

import os


class ArtIllumina(Source):
	"""Class for the ART Illumina (http://www.niehs.nih.gov/research/resources/software/biostatistics/art/).

	Single-end reads and pair-end reads simulations are supported. For pair-end simulations,
	lengths of both ends must be equal.

	Args:
		fasta (str): File name of the genome from which read tuples are created (FASTA file). Corresponding ART parameter: ``-i --in``.
		sequences (set of int or str): FASTA sequences to extract. Sequences can be specified either by their ids, or by their names.
		coverage (float): Average coverage of the genome. Corresponding ART parameter: ``-f --fcov``.
		number_of_read_tuples (int): Number of read tuples.
		read_length_1 (int): Length of the first read. Corresponding ART parameter: ``-l --len``.
		read_length_2 (int): Length of the second read (if zero, then single-end reads are simulated). Corresponding ART parameter: ``-l --len``.
		distance (int): Mean inner distance between reads. Corresponding ART parameter: ``-m --mflen``.
		distance_deviation (int): Standard devation of inner distances between reads. Corresponding ART parameter: ``-s --sdev``.
		rng_seed (int): Seed for simulator's random number generator. Corresponding ART parameter: ``-rs --rndSeed``.
		other_params (str): Other parameters which are used on commandline.

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
			other_params="",
	):

		if read_length_2 == 0:
			ends = 1
		else:
			ends = 2
			self.distance = distance
			self.distance_deviation = distance_deviation
			if read_length_1 != read_length_2:
				rnftools.utils.error(
					"art_illumina can simulate only pairs with equal lengths",
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
				exception=ValueError,
			)

		self.number_of_read_tuples = number_of_read_tuples
		self.coverage = coverage

		self.art_prefix = os.path.join(
			self.get_dir(),
			"art_illumina_files.{}.{}".format("se" if self.number_of_read_tuples == 1 else "pe", self.genome_id)
		)

		self._sam1_fn = self.art_prefix + ".sam"
		# self._sam2_fn = self.art_prefix+".sam"
		self._sam2_fn = self.art_prefix + ".corrected.sam"

	def get_input(self):
		return [
			self._fa_fn,
			self._fai_fn,
		]

	def get_output(self):
		if self._reads_in_tuple == 1:
			fqs = [self.art_prefix + ".fq"]
		else:
			fqs = [self.art_prefix + "1.fq", self.art_prefix + "2.fq"]
		return [
				   self._fq_fn,
				   self._sam1_fn,
				   self._sam2_fn,
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
				paired_params = "-p -m {dist} -s {dist_dev}".format(
					dist=self.distance,
					dist_dev=self.distance_deviation,
				)
			else:
				paired_params = ""

			command_1 = """
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
				art_il="art_illumina",
				paired_params=paired_params,
				fasta=self._fa_fn,
				rlen=self.read_length_1,
				other_params=self.other_params,
				coverage=self.coverage,
				o_pref=self.art_prefix,
				rng_seed=self._rng_seed,
			)

			# correction of header (bug in ART)
			command_2 = """
				cat "{sam_1}" | \
				grep -v ^@ | \
				"{samtools}" view -h -T "{fa}" - \
				> "{sam_2}"
			""".format(
				samtools="samtools",
				sam_1=self._sam1_fn,
				sam_2=self._sam2_fn,
				fa=self._fa_fn,
			)

			rnftools.utils.shell(command_1)
			rnftools.utils.shell(command_2)

			with open(self._fq_fn, "w+") as fq_fo:
				with open(self._fai_fn) as fai_fo:
					self.recode_sam_reads(
						sam_fn=self._sam2_fn,
						fastq_rnf_fo=fq_fo,
						fai_fo=fai_fo,
						genome_id=self.genome_id,
						number_of_read_tuples=10 ** 9,
						simulator_name="art-illumina",
						allow_unmapped=False,
					)

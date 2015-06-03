import snakemake
import smbl

import abc
import re
import os
import pysam

import rnftools

class Source(object):
	"""	Abstract class for a genome from which read tuples are simulated.

	Args:
		fasta (str): File name of the genome from which reads are created (FASTA file).
		reads_in_tuple (int): Number of reads in each read tuple.
		rng_seed (int): Seed for simulator's random number generator.
		number_of_required_cores (int): Number of cores used by the simulator. It can be used to prevent running other threads at the same time.
	"""

	__metaclass__ = abc.ABCMeta

	def __init__(self,
				fasta,
				reads_in_tuple,
				rng_seed,
				number_of_required_cores=1
			):
		rnftools.mishmash.add_source(self)
		self._rng_seed = rng_seed
		self._reads_in_tuple=reads_in_tuple

		self._sample=rnftools.mishmash.current_sample()
		self._sample.add_source(self)
		self.genome_id=len(self._sample.get_sources())
		self._number_of_required_cores=number_of_required_cores

		self._name=str(self.genome_id).zfill(3)

		self._dir=os.path.join(
					self._sample.get_dir(),
					self._name
				)
		self._fa_fn=fasta
		self._fai_fn = fasta+".fai"
		self._fq_fn=os.path.join(self._dir,"_final_reads.fq")
		self.dict_chr_ids = {}
		self.dict_chr_lengths = {}

	############################################################################
	############################################################################

	def get_dir(self):
		"""Get working directory.

		Returns:
			str: Working directory.
		"""
		return self._dir

	def get_genome_id(self):
		"""Get genome ID.

		Returns:
			int: Genome ID.
		"""
		return self.genome_id

	def get_reads_in_tuple(self):
		"""Get number of entries in a read tuple.

		Returns:
			int: Number of reads in a read tuple.
		"""
		return self._reads_in_tuple

	def get_number_of_required_cores(self):
		"""Get number of required cores.

		Returns:
			int: Number of required cores.
		"""
		return self._number_of_required_cores

	def clean(self):
		"""Clean working directory.
		"""
		snakemake.shell('rm -fR "{}"'.format(self.get_dir()))

	############################################################################
	############################################################################

	def fa_fn(self):
		"""Get input FASTA file. It can be an empty list.

		Returns:
			str: Input FASTA file.
		"""
		return self._fa_fn

	def fq_fn(self):
		"""Get file name of the output FASTQ file.

		Returns:
			str: Output FASTQ file
		"""
		return self._fq_fn


	@abc.abstractmethod
	def get_input(self):
		"""Get list of input files (required to do simulation).

		Returns:
			list: List of input files
		"""
		return

	############################################################################
	############################################################################

	@abc.abstractmethod
	def get_output(self):
		"""Get list of output files (created during simulation).
		
		Returns:
			list: List of input files
		"""
		return

	############################################################################
	############################################################################

	@abc.abstractmethod
	def create_fq(self):
		"""Simulate reads.
		"""
		return

	# todo: make this method static
	# todo: check if it can work with a bam file
	# todo: can it work as a pipe?
	@staticmethod
	def recode_sam_reads(
				sam,
				fastq,
				fai,
				genome_id,
				number_of_read_tuples=10**9,
				simulator_name=None,
				allow_unmapped=False,
			):
		"""Create FASTQ file from SAM file.

		Args:
			sam (str): Name of SAM file.
			number_of_read_tuples (int): Expected number of read tuples (to set width of read tuple id).
			simulator_name (str): Name of the simulator. Used for comment in read tuple name.

		Raises:
			NotImplementedError
		"""

		fai_index = FaiIndex(fai)
		#last_read_tuple_name=[]
		read_tuple_id_width=len(format(number_of_read_tuples,'x'))
		fq_creator=rnftools.rnfformat.FqCreator(
					fastq=fastq,
					read_tuple_id_width=read_tuple_id_width,
					genome_id_width=2,
					chr_id_width=fai_index.chr_id_width,
					coor_width=fai_index.coor_width,
					info_reads_in_tuple=True,
					info_simulator=simulator_name,					
				)

		#todo: check if clipping corrections is well implemented
		cigar_reg_shift=re.compile("([0-9]+)([MDNP=X])")

		#todo: other upac codes
		reverse_complement_dict = {
			"A":"T",
			"T":"A",
			"C":"G",
			"G":"C",
			"N":"N",
		}

		read_tuple_id=0
		last_read_tuple_name=None
		with pysam.AlignmentFile(sam, "rb") as samfile:
			for alignment in samfile:
				if alignment.query_name!=last_read_tuple_name and last_read_tuple_name is not None:
					read_tuple_id+=1
				last_read_tuple_name = alignment.query_name

				if alignment.is_unmapped:
					smbl.messages.error(
						"SAM files used for conversion should not contain unaligned segments. "
						"This condition is broken by read tuple "
						"'{}' in file '{}'.".format(alignment.query_name,sam),
						program="RNFtools",
						subprogram="MIShmash",
						exception=NotImplementedError
					)


				if alignment.is_reverse:
					direction  = "R"
					bases      = "".join([reverse_complement_dict[nucl]
						for nucl in alignment.seq[::-1]])
					qualities  = alignment.qual[::-1]
				else:
					direction  = "F"
					bases      = alignment.seq[:]
					qualities  = alignment.qual[:]

				# todo: are chromosomes in bam sorted correctly (the same order as in FASTA)?
				if fai_index.dict_chr_ids!={}:
					chr_id=fai_index.dict_chr_ids[ samfile.getrname(alignment.reference_id) ]
				else:
					chr_id="0"

				left=int(alignment.reference_start)+1
				right=left-1
				for (steps,operation) in cigar_reg_shift.findall(alignment.cigarstring):
					right+=int(steps)

				segment=rnftools.rnfformat.Segment(
						genome_id=genome_id,
						chr_id=chr_id,
						direction=direction,
						left=left,
						right=right,
					)

				fq_creator.add_read(
						read_tuple_id=read_tuple_id,
						bases=bases,
						qualities=qualities,
						segments=[segment],
					)
		fq_creator.flush_read_tuple()

	"""Load dictionaries with sizes of chromosomes and with id-name correspondance.
	"""
	def load_fai(self):
		self.dict_chr_ids = {}
		self.dict_chr_lengths = {}

		# parsing FAI file
		with open(self._fai_fn) as f:
			"""
			   FAI format

			1) the name of the sequence
			2) the length of the sequence
			3) the offset of the first base in the file
			4) the number of bases in each fasta line
			5) the number of bytes in each fasta line
			"""

			i=1
			for line in f:
				if line.strip()!="":
					parts=line.split("\t")
					chr=parts[0]
					chr_len=int(parts[1])
					self.dict_chr_ids[chr]=i
					self.dict_chr_lengths[chr]=chr_len
					i+=1

		self.number_of_chromosomes=len(self.dict_chr_ids)
		self.chr_id_width=len(str(self.number_of_chromosomes))
		self.coor_width=len(str(max(self.dict_chr_lengths.values())))

class FaiIndex:
	def __init__(self, fai):
		self.dict_chr_ids = {}
		self.dict_chr_lengths = {}

		# parsing FAI file
		with open(fai) as f:
			"""
			   FAI format

			1) the name of the sequence
			2) the length of the sequence
			3) the offset of the first base in the file
			4) the number of bases in each fasta line
			5) the number of bytes in each fasta line
			"""

			i=1
			for line in f:
				if line.strip()!="":
					parts=line.split("\t")
					chr=parts[0]
					chr_len=int(parts[1])
					self.dict_chr_ids[chr]=i
					self.dict_chr_lengths[chr]=chr_len
					i+=1

		self.number_of_chromosomes=len(self.dict_chr_ids)
		self.chr_id_width=len(str(self.number_of_chromosomes))
		self.coor_width=len(str(max(self.dict_chr_lengths.values())))



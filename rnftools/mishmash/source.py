import snakemake
import smbl

import abc
import re
import os
import pysam

import rnftools

class Source(object):
	"""	Abstract class for a genome from which reads are simulated.

	Args:
		fasta (str): File name of the genome from which reads are created (FASTA file).
		reads_in_tuple (int): Number of reads in each read tuple.
		rng_seed (int): Seed for simulator's random number generator.
		number_of_required_cores (int): Number of cores used by the simulator. It can be used to prevent running other threads at the same time.
	"""

	__metaclass__ = abc.ABCMeta

	def __init__(self, fasta, reads_in_tuple, rng_seed, number_of_required_cores=1):
		rnftools.mishmash.add_source(self)
		self._rng_seed = rng_seed
		self._reads_in_tuple=reads_in_tuple

		self._sample=rnftools.mishmash.current_sample()
		self._sample.add_source(self)
		self.source_id=len(self._sample.get_sources())
		self._number_of_required_cores=number_of_required_cores

		self._name=str(self.source_id).zfill(3)

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

	def get_source_id(self):
		"""Get genome ID.

		Returns:
			int: Genome ID.
		"""
		return self.source_id

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

	def _fq_buffer(self,read_id,segments_buffer,sequences_buffer,rn_formatter,simulator_name=""):
		"""From local buffers, create FASTQ string.

		Args:
			reads_id (int): ID of read tuple.
			segments_buffer(list of rnftools.rnfformat.Segment): Buffer of segments.
			sequences_buffer(list of (bases,qualities)):  Buffer of sequences.
			rn_formatter (rnftools.rnfformat.RnFormatter): Read formatter.
			simulator_name (str): Name of the simulator. Used for comment in read name.

		Returns:
			str: Part of FASTQ file.
		"""
		read_suffix_comment_buffer=[]
		if len(segments_buffer)==1:
			read_suffix_comment_buffer.append("single-end")
		elif len(segments_buffer)==2 and set([segments_buffer[i].direction for i in [0,1]])==set(["R","F"]):
			read_suffix_comment_buffer.append("paired-end")

		if simulator_name!="":
			read_suffix_comment_buffer.append(simulator_name)

		if len(read_suffix_comment_buffer)!=0:
			read_suffix="[{}]".format(",".join(read_suffix_comment_buffer))
		else:
			read_suffix=""

		rnf_read = rnftools.rnfformat.Read(segments=segments_buffer,read_id=read_id,suffix=read_suffix)
		rnf_read_name = rn_formatter.process_read(read=rnf_read)
		to_return = [
			"".
				join([
					"@",rnf_read_name,"/{}".format(str(b_i)) if len(sequences_buffer)>1 else "",
						os.linesep,
					sequences_buffer[b_i-1][0],
						os.linesep,
					"+",
						os.linesep,
					sequences_buffer[b_i-1][1],
						os.linesep,
				])
			for b_i in range(1,len(sequences_buffer)+1)
		]
		return "".join(to_return)

	def recode_sam_reads(self,sam,number_of_read_tuples=10**9,simulator_name=""):
		"""Create FASTQ file from SAM file.

		Args:
			sam (str): Name of SAM file.
			number_of_read_tuples (int): Number of read tuples. It is needed to set width of read tuple id.
			simulator_name (str): Name of the simulator. Used for comment in read name.

		Raises:
			NotImplementedError
		"""
		self.load_fai()
		last_read_name=[]
		read_id=0

		id_str_size=len(format(number_of_read_tuples,'x'))

		rn_formatter = rnftools.rnfformat.RnFormatter(
				id_str_size=id_str_size,
				source_str_size=2,
				chr_str_size=self.chr_str_size,
				pos_str_size=self.pos_str_size
			)

		segments_buffer=[]
		sequences_buffer=[]
		last_read_name=""

		reverse_complement_dict = {
			"A":"T",
			"T":"A",
			"C":"G",
			"G":"C",
			"N":"N"
		}

		cigar_reg_shift=re.compile("([0-9]+)([MDNP=X])")
		with pysam.AlignmentFile(sam, "rb") as samfile:
			with open(self._fq_fn, "w+") as fqfile:
				for alignment in samfile:
					if alignment.query_name!=last_read_name and last_read_name!="":
						read_id+=1
						fqfile.write(
							self._fq_buffer(
									read_id=read_id,
									segments_buffer=segments_buffer,
									sequences_buffer=sequences_buffer,
									rn_formatter=rn_formatter,
									simulator_name=simulator_name,
								)
						)
						segments_buffer = []
						sequences_buffer = []
					last_read_name = alignment.query_name

					if alignment.is_unmapped:
						smbl.messages.error(
							"SAM files used for conversion should not contain unaligned segments. "
							"This condition is broken by read "
							"'{}' in file '{}'.".format(alignment.query_name,sam),program="RNFtools",subprogram="MIShmash",exception=NotImplementedError)


					if alignment.is_reverse:
						direction  = "R"
						bases      = "".join([reverse_complement_dict[nucl]
							for nucl in alignment.seq[::-1]])
						qualities  = alignment.qual[::-1]
					else:
						direction  = "F"
						bases      = alignment.seq[:]
						qualities  = alignment.qual[:]

					if self.dict_chr_ids!={}:
						chr_id=self.dict_chr_ids[ samfile.getrname(alignment.reference_id) ]
					else:
						chr_id="0"

					left=int(alignment.reference_start)+1
					right=left-1
					for (steps,operation) in cigar_reg_shift.findall(alignment.cigarstring):
						right+=int(steps)

					segment=rnftools.rnfformat.Segment(
							source=self.source_id,
							chr=chr_id,
							direction=direction,
							left=left,
							right=right
						)

					segments_buffer.append(segment)
					sequences_buffer.append( (bases,qualities) )

				fqfile.write(
					self._fq_buffer(
							read_id=read_id,
							segments_buffer=segments_buffer,
							sequences_buffer=sequences_buffer,
							rn_formatter=rn_formatter,
							simulator_name=simulator_name,
						)
				)


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
		self.chr_str_size=len(str(self.number_of_chromosomes))
		self.pos_str_size=len(str(max(self.dict_chr_lengths.values())))

import snakemake
import smbl
import rnftools
import os

class FqCreator:
	"""Class for writing RNF reads to FASTQ files.

	Args:
		fastq_fo (str): Output FASTQ file - file object.
		read_tuple_id_width (int): Maximal expected string length of read tuple ID.
		genome_id_width (int): Maximal expected string length of genome ID.
		chr_id_width (int): Maximal expected string length of chromosome ID.
		coor_width (int): Maximal expected string length of a coordinate.
		info_reads_in_tuple (bool): Include information about reads as a RNF comment.
		info_simulator (str): Name of used simulator (to be included as a RNF comment).
	"""

	def __init__(
				self,
				fastq_fo,
				read_tuple_id_width=16,
				genome_id_width=2,
				chr_id_width=2,
				coor_width=8,
				info_reads_in_tuple=True,
				info_simulator=None,
			):
		self._info_simulator=info_simulator
		self._info_reads_in_tuple=info_reads_in_tuple
		self._rnf_profile=rnftools.rnfformat.RnfProfile(
				read_tuple_id_width=read_tuple_id_width,
				genome_id_width=genome_id_width,
				chr_id_width=chr_id_width,
				coor_width=coor_width,
			)
		self._fq_file=fastq_fo
		self.current_read_tuple_id=None
		self.empty()

	def flush_read_tuple(self):
		"""Flush the internal buffer of reads.
		"""
		#print("flush")
		suffix_comment_buffer=[]
		if self._info_simulator is not None:
			suffix_comment_buffer.append(self._info_simulator)
		if self._info_reads_in_tuple:
			#todo: orientation (FF, FR, etc.)
			#orientation="".join([])
			suffix_comment_buffer.append("reads-in-tuple:{}".format(len(self.seqs_bases)))
		if len(suffix_comment_buffer)!=0:
			suffix_comment="[{}]".format(",".join(suffix_comment_buffer))
		else:
			suffix_comment=""

		rnf_name=self._rnf_profile.get_rnf_name(
					rnftools.rnfformat.ReadTuple(
							segments=self.segments,
							read_tuple_id=self.current_read_tuple_id,
							suffix=suffix_comment,
						)
				)
		fq_reads = [
			os.linesep.join([
					"@{rnf_name}{read_suffix}".format(
							rnf_name=rnf_name,
							read_suffix="/{}".format(str(i+1)) if len(self.seqs_bases)>1 else "",
						),
					self.seqs_bases[i],
					"+",
					self.seqs_qualities[i],
				])
			for i in range(len(self.seqs_bases))
		]
		self._fq_file.write(os.linesep.join(fq_reads))
		self._fq_file.write(os.linesep)
		self.empty()

	def empty(self):
		"""Empty all internal buffers.
		"""
		self.seqs_bases=[]
		self.seqs_qualities=[]
		self.segments=[]

	def add_read(self,
				read_tuple_id,
				bases,
				qualities,
				segments,
			):

		"""Add a new read to the current buffer. If it is a new read tuple, the buffer will be flushed.

		Args:
			read_tuple_id (int): ID of the read tuple.
			bases (str): Sequence of bases.
			qualities (str): Sequence of FASTQ qualities.
			segments (list of rnftools.rnfformat.segment): List of segments constituting the read.
		"""

		assert type(bases) is str, "Wrong type of bases: '{}'".format(bases)
		assert type(qualities) is str, "Wrong type of qualities: '{}'".format(qualities)
		assert type(segments) is tuple or type(segments) is list

		if self.current_read_tuple_id!=read_tuple_id and self.current_read_tuple_id is not None:
			self.flush_read_tuple()
		self.current_read_tuple_id=read_tuple_id

		self.seqs_bases.append(bases)
		self.seqs_qualities.append(qualities)
		self.segments.extend(segments)

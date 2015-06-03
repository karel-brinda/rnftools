import snakemake
import smbl
import rnftools
import os

class FqCreator:
	"""Class for writing RNF reads to FASTQ files.

	Args:
		fastq (str): Output FASTQ file.
		read_tuple_id_width (int): Maximal expected string length of read tuple ID.
		genome_id_width (int): Maximal expected string length of genome ID.
		chr_id_width (int): Maximal expected string length of chromosome ID.
		coor_width (int): Maximal expected string length of a coordinate.
		info_reads_in_tuple (bool): Include information about reads as a RNF comment.
		info_simulator (str): Name of used simulator (to be included as a RNF comment).
	"""

	def __init__(
				self,
				fastq,
				read_tuple_id_width=16,
				genome_id_width=2,
				chr_id_width=2,
				coor_width=8,
				info_reads_in_tuple=True,
				info_simulator=None,
			):
		self._fq=fastq
		self._info_simulator=info_simulator
		self._info_reads_in_tuple=info_reads_in_tuple
		self._formatter=rnftools.rnfformat.RnFormatter(
				read_tuple_id_width=read_tuple_id_width,
				genome_id_width=genome_id_width,
				chr_id_width=chr_id_width,
				coor_width=coor_width,
			)
		self._fq_file=open(self._fq,"w+")
		self.current_read_tuple_id=None
		self.empty()

	def __del__(self):
		self.flush_read_tuple()
		self._fq_file.close()

	def flush_read_tuple(self):
		"""Flush the internal buffer of reads.
		"""
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

		rnf_name=self._formatter.get_rnf_name(
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

#
#	def _fq_buffer(self,read_tuple_id,segments_buffer,sequences_buffer,rn_formatter,simulator_name=""):
#		"""From local buffers, create FASTQ string.
#
#		Args:
#			read_tuple_id (int): ID of read tuple.
#			segments_buffer(list of rnftools.rnfformat.Segment): Buffer of segments.
#			sequences_buffer(list of (bases,qualities)):  Buffer of sequences.
#			rn_formatter (rnftools.rnfformat.RnFormatter): Read formatter.
#			simulator_name (str): Name of the simulator. Used for comment in tuple read name.
#
#		Returns:
#			str: Part of FASTQ file.
#		"""
#		read_tuple_suffix_comment_buffer=[]
#		if len(segments_buffer)==1:
#			read_tuple_suffix_comment_buffer.append("single-end")
#		elif len(segments_buffer)==2 and set([segments_buffer[i].direction for i in [0,1]])==set(["R","F"]):
#			read_tuple_suffix_comment_buffer.append("paired-end")
#
#		if simulator_name!="":
#			read_tuple_suffix_comment_buffer.append(simulator_name)
#
#		if len(read_tuple_suffix_comment_buffer)!=0:
#			read_tuple_suffix="[{}]".format(",".join(read_tuple_suffix_comment_buffer))
#		else:
#			read_tuple_suffix=""
#
#		rnf_read_tuple = rnftools.rnfformat.ReadTuple(segments=segments_buffer,read_tuple_id=read_tuple_id,suffix=read_tuple_suffix)
#		rnf_read_tuple_name = rn_formatter.process_read_tuple(read_tuple=rnf_read_tuple)
#		to_return = [
#			"".
#				join([
#					"@",rnf_read_tuple_name,"/{}".format(str(b_i)) if len(sequences_buffer)>1 else "",
#						os.linesep,
#					sequences_buffer[b_i-1][0],
#						os.linesep,
#					"+",
#						os.linesep,
#					sequences_buffer[b_i-1][1],
#						os.linesep,
#				])
#			for b_i in range(1,len(sequences_buffer)+1)
#		]
#		return "".join(to_return)
class RnFormatter:
	"""Class for formatting RNF reads in the same format.

	Args:
		read_tuple_id_width (int): Maximal expected string length of read tuple ID.
		genome_id_width (int): Maximal expected string length of genome ID.
		chr_id_width (int): Maximal expected string length of chromosome ID.
		coor_width (int): Maximal expected string length of a coordinate.
	"""

	def __init__(
				self,
				read_tuple_id_width=16,
				genome_id_width=2,
				chr_id_width=2,
				coor_width=8,
			):

		self.read_tuple_id_width=read_tuple_id_width
		self.genome_id_width=genome_id_width
		self.chr_id_width=chr_id_width
		self.coor_width=coor_width
	
	def process_read_tuple(self,read_tuple):
		"""Get well-formatted RNF representation of a read.

		read_tuple (rnftools.rnFormatter.ReadTuple): Read tuple.
		"""
		return read_tuple.stringize(
					read_tuple_id_width=self.read_tuple_id_width,
					genome_id_width=self.genome_id_width,
					chr_id_width=self.chr_id_width,
					coor_width=self.coor_width
				)

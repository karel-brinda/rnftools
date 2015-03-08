class RnFormatter:
	"""Class for formatting RNF reads in the same format.
	"""

	def __init__(
				self,
				id_str_size=16,
				source_str_size=2,
				chr_str_size=2,
				pos_str_size=8,
			):
		"""
		:param id_str_size: Maximal expected length of string of ID of a read.
		:type  id_str_size: int
		:param source_str_size: Maximal expected length of string of ID of genome.
		:type  source_str_size: int
		:param chr_str_size: Maximal expected length of string of ID of chromosome.
		:type  chr_str_size: int
		:param pos_str_size: Maximal expected length of string of maximal coordinate.
		:type  pos_str_size: int
		"""

		self.id_str_size=id_str_size
		self.source_str_size=source_str_size
		self.chr_str_size=chr_str_size
		self.pos_str_size=pos_str_size
	
	def process_read(self,read):
		"""Get well-formatted RNF representation of a read.

		:param read: Read.
		:type  read: rnftools.rnFormatter.Read
		"""
		return read.stringize(
					id_str_size=self.id_str_size,
					source_str_size=self.source_str_size,
					chr_str_size=self.chr_str_size,
					pos_str_size=self.pos_str_size
				)

class RnFormatter:
	"""Class for formatting RNF reads in the same format.

	Args:
		id_str_size (int): Maximal expected length of string of ID of a read.
		source_str_size (int): Maximal expected length of string of ID of genome.
		chr_str_size (int): Maximal expected length of string of ID of chromosome.
		pos_str_size (int): Maximal expected length of string of maximal coordinate.
	"""

	def __init__(
				self,
				id_str_size=16,
				source_str_size=2,
				chr_str_size=2,
				pos_str_size=8,
			):

		self.id_str_size=id_str_size
		self.source_str_size=source_str_size
		self.chr_str_size=chr_str_size
		self.pos_str_size=pos_str_size
	
	def process_read(self,read):
		"""Get well-formatted RNF representation of a read.

		read (rnftools.rnFormatter.Read): Read.
		"""
		return read.stringize(
					id_str_size=self.id_str_size,
					source_str_size=self.source_str_size,
					chr_str_size=self.chr_str_size,
					pos_str_size=self.pos_str_size
				)

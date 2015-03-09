import re

segment_destr_pattern = re.compile(r'^\(([0-9]+),([0-9]+),([FRN]),([0-9]+),([0-9]+)\)$')

class Segment:
	"""Class for a single segment in a RNF read name.
	"""

	def __init__(
				self,
				source=0,
				chr=0,
				direction="N",
				left=0,
				right=0
			):
		"""
		:param source: ID of source.
		:type  fa: int
		:param chr: ID of chromosome.
		:type  chr: int
		:param direction: Direction (F/R/N).
		:type  direction: str
		:param left: Leftmost coordinate.
		:type  left: int
		:param right: Rightmost coordinate.
		:type  right: int
		:raises: ValueError
		"""

		self.source=int(source)
		self.chr=int(chr)
		self.direction=direction
		self.left=int(left)
		self.right=int(right)

		if not (self.right == 0 or self.right >= self.left):
			raise valueError("Leftmost coordinate cannot be higher than rightmost coordinate")

	def stringize(
				self,
				source_str_size=1,
				chr_str_size=1,
				pos_str_size=1
			):
		"""Create RNF representation of this segment.

		:param source_str_size: Maximal expected length of string of ID of genome.
		:type  source_str_size: int
		:param chr_str_size: Maximal expected length of string of ID of chromosome.
		:type  chr_str_size: int
		:param pos_str_size: Maximal expected length of string of maximal coordinate.
		:type  pos_str_size: int
		"""

		pos_str_size=max(pos_str_size,len(str(self.left)),len(str(self.right)))
		return "({},{},{},{},{})".format(
				str(self.source).zfill(source_str_size),
				str(self.chr).zfill(chr_str_size),
				self.direction,
				str(self.left).zfill(pos_str_size),
				str(self.right).zfill(pos_str_size)
			)

	def destringize(self,string):
		"""Get RNF values for this segment from its textual representation and
		save them into this object.

		:param string: Textual representation of a segment.
		:type  string: str
		"""

		m=segment_destr_pattern.match(string)
		self.source=int(m.group(1))
		self.chr=int(m.group(2))
		self.direction=m.group(3)
		self.left=int(m.group(4))
		self.right=int(m.group(5))
		#TODO: checks:
		#   1. if reg match ok
		#   2. if variables ok left <= right, etc.

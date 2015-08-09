import re
from .RnfProfile import RnfProfile

segment_destr_pattern = re.compile(r'^\(([0-9]+),([0-9]+),([FRN]),([0-9]+),([0-9]+)\)$')

class Segment:
	"""Class for a single segment in a RNF read name.
	"""

	def __init__(
				self,
				genome_id=0,
				chr_id=0,
				direction="N",
				left=0,
				right=0
			):
		"""
		Args:
			genome_id (int): ID of genome.
			chr_id (int): ID of chromosome.
			direction (str): Direction (F/R/N).
			left (int): Leftmost coordinate.
			right (int): Rightmost coordinate.
		
		Raises:
			ValueError
		"""

		self.genome_id=int(genome_id)
		self.chr_id=int(chr_id)
		self.direction=direction
		self.left=int(left)
		self.right=int(right)

		if not (self.right == 0 or self.right >= self.left):
			smbl.messages.error("Leftmost coordinate cannot be higher than rightmost coordinate.",program="RNFtools",subprogram="RNF format",exception=ValueError)

	def stringize(
				self,
				rnf_profile=RnfProfile(),
			):
		"""Create RNF representation of this segment.

		Args:
			rnf_profile (rnftools.rnfformat.RnfProfile): RNF profile (with widths).
		"""

		coor_width=max(rnf_profile.coor_width,len(str(self.left)),len(str(self.right)))
		return "({},{},{},{},{})".format(
				str(self.genome_id).zfill(rnf_profile.genome_id_width),
				str(self.chr_id).zfill(rnf_profile.chr_id_width),
				self.direction,
				str(self.left).zfill(coor_width),
				str(self.right).zfill(coor_width)
			)

	def destringize(self,string):
		"""Get RNF values for this segment from its textual representation and
		save them into this object.

		Args:
			string (str): Textual representation of a segment.
		"""

		m=segment_destr_pattern.match(string)
		self.genome_id=int(m.group(1))
		self.chr_id=int(m.group(2))
		self.direction=m.group(3)
		self.left=int(m.group(4))
		self.right=int(m.group(5))

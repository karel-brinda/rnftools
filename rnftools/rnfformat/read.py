import rnftools.rnfformat
import re

read_destr_pattern = re.compile(r'(.*)__([0-9abcdef]+)__(\([0-9abcdefFRN,]*\))(,\([0-9abcdefFRN,]*\))*__(.*)')

class Read:
	"""Class for a RNF read.
	"""

	def __init__(
				self,
				segments=[],
				read_id=0,
				prefix="",
				suffix="",
			):
		"""
		:param segments: Segments of the read.
		:type  segments: list of rnftools.rnfformat.Segment
		:param read_id: Read ID.
		:type  read_id: int
		:param prefix: Prefix for the read name.
		:type  prefix: str
		:param suffix: Suffix for the read name.
		:type  suffix: str
		"""
		self.read_id=read_id
		self.segments=segments
		self.prefix=prefix
		self.suffix=suffix

	def stringize(
				self,
				id_str_size=1,
				source_str_size=1,
				chr_str_size=1,
				pos_str_size=1
			):
		"""Create RNF representation of this read.

		:param id_str_size: Maximal expected length of string of ID of a read.
		:type  id_str_size: int
		:param source_str_size: Maximal expected length of string of ID of genome.
		:type  source_str_size: int
		:param chr_str_size: Maximal expected length of string of ID of chromosome.
		:type  chr_str_size: int
		:param pos_str_size: Maximal expected length of string of maximal coordinate.
		:type  pos_str_size: int
		"""

		sorted_segments = sorted(self.segments,
								key=lambda x: (
									x.source*(10**23) + 
									x.chr*(10**21) +
									(x.left+(int(x.left==0)*x.right-1))*(10**11) +
									x.right*(10**1) + 
									int(x.direction=="F")
								)
							)

		segments_strings=[
				x.stringize(
					source_str_size=source_str_size,
					chr_str_size=chr_str_size,
					pos_str_size=pos_str_size
				) for x in sorted_segments
			]

		read_name="__".join(
			[
				self.prefix,
				format(self.read_id,'x').zfill(id_str_size),
				",".join(segments_strings),
				self.suffix,
			]
		)

		return read_name

	def destringize(self, string):
		"""Get RNF values for this read from its textual representation and save them 
		into this object.

		:param string: Textual representation of a read.
		:type  string: str
		"""

		#todo: assert -- starting with (, ending with )
		#(prefix,read_id,segments_t,suffix)=(text).split("__")
		#segments=segments_t.split("),(")
		m=read_destr_pattern.match(string)
		if not m:
			raise ValueError("'{}' is not a valid read name with respect to the RNF specification".format(string))
		groups=m.groups()
		#todo: check number of groups
		self.prefix=groups[0]
		read_id=groups[1]
		self.read_id=int(read_id,16)
		self.segments=[]
		segments_str=groups[2:-1]
		for b_str in segments_str:
			if b_str is not None:
				if b_str[0]==",":
					b_str=b_str[1:]
				b=rnftools.rnfformat.Segment()
				b.destringize(b_str)
				self.segments.append(b)
		self.suffix=groups[-1]

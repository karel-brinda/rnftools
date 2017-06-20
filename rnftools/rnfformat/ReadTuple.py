import rnftools.rnfformat
from .RnfProfile import RnfProfile
import re

read_tuple_destr_pattern = re.compile(r'(.*)__([0-9abcdef]+)__(\([0-9abcdefFRN,]*\))(,\([0-9abcdefFRN,]*\))*__(.*)')


class ReadTuple:
	"""Class for a RNF read tuple.

	Args:
		segments (list of rnftools.rnfformat.Segment): Segments of the read.
		read_tuple_id (int): Read tuple ID.
		prefix (str): Prefix for the read name.
		suffix (str): Suffix for the read name.
	"""

	def __init__(
			self,
			segments=[],
			read_tuple_id=0,
			prefix="",
			suffix="",
	):

		assert type(segments) is tuple or type(segments) is list, "Wrong type of segments: '{}'".format(segments)
		assert type(read_tuple_id) is int, "Wrong type of read_tuple_id: '{}'".format(read_tuple_id)
		assert type(prefix) is str, "Wrong type of prefix: '{}'".format(prefix)
		assert type(suffix) is str, "Wrong type of suffix: '{}'".format(suffix)

		self.read_tuple_id = read_tuple_id
		self.segments = segments
		self.prefix = prefix
		self.suffix = suffix

	def stringize(
			self,
			rnf_profile=RnfProfile(),
	):
		"""Create RNF representation of this read.

		Args:
			read_tuple_id_width (int): Maximal expected string length of read tuple ID.
			genome_id_width (int): Maximal expected string length of genome ID.
			chr_id_width (int): Maximal expected string length of chromosome ID.
			coor_width (int): Maximal expected string length of a coordinate.
		"""

		sorted_segments = sorted(self.segments,
			key=lambda x: (
				x.genome_id * (10 ** 23) +
				x.chr_id * (10 ** 21) +
				(x.left + (int(x.left == 0) * x.right - 1)) * (10 ** 11) +
				x.right * (10 ** 1) +
				int(x.direction == "F")
			)
		)

		segments_strings = [
			x.stringize(rnf_profile) for x in sorted_segments
		]

		read_tuple_name = "__".join(
			[
				self.prefix,
				format(self.read_tuple_id, 'x').zfill(rnf_profile.read_tuple_id_width),
				",".join(segments_strings),
				self.suffix,
			]
		)

		return read_tuple_name

	def destringize(self, string):
		"""Get RNF values for this read from its textual representation and save them 
		into this object.

		Args:
			string(str): Textual representation of a read.

		Raises:
			ValueError
		"""

		# todo: assert -- starting with (, ending with )
		# (prefix,read_tuple_id,segments_t,suffix)=(text).split("__")
		# segments=segments_t.split("),(")
		m = read_tuple_destr_pattern.match(string)
		if not m:
			smbl.messages.error("'{}' is not a valid read name with respect to the RNF specification".format(string),
				program="RNFtools", subprogram="RNF format", exception=ValueError)
		groups = m.groups()
		# todo: check number of groups
		self.prefix = groups[0]
		read_tuple_id = groups[1]
		self.read_tuple_id = int(read_tuple_id, 16)
		self.segments = []
		segments_str = groups[2:-1]
		for b_str in segments_str:
			if b_str is not None:
				if b_str[0] == ",":
					b_str = b_str[1:]
				b = rnftools.rnfformat.Segment()
				b.destringize(b_str)
				self.segments.append(b)
		self.suffix = groups[-1]

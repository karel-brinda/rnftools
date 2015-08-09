class RnfProfile:
	"""Class for profile of RNF reads (widths).

	Args:
		prefix_width (int): Length of prefix.
		read_tuple_id_width (int): Width of read tuple ID
		genome_id_width (int): Width of genome ID.
		chr_id_width (int): Width of chromosome ID.
		coor_width (int): Width of coordinate width.
		read_tuple_name (str): Read tuple name to initialize all the values.

	Attributes:
		prefix_width (int): Length of prefix.
		read_tuple_id_width (int): Width of read tuple ID
		genome_id_width (int): Width of genome ID.
		chr_id_width (int): Width of chromosome ID.
		coor_width (int): Width of coordinate width.
	"""

	def __init__(self,
				prefix_width=0,
				read_tuple_id_width=8,
				genome_id_width=1,
				chr_id_width=2,
				coor_width=9,
				read_tuple_name=None
			):
		self.prefix_width=prefix_width
		self.read_tuple_id_width=read_tuple_id_width
		self.genome_id_width=genome_id_width
		self.chr_id_width=chr_id_width
		self.coor_width=coor_width

		if read_tuple_name is not None:
			self.load(read_tuple_name)

	def __str__(self):
		return str(list([
				self.prefix_width,
				self.read_tuple_id_width,
				self.genome_id_width,
				self.chr_id_width,
				self.coor_width,
			]))
		#return "{{prefix_width:{},read_tuple_id_width:{},genome_id_width:{},chr_id_width:{},coor_width:{}}}".format(
		#		self.prefix_width,
		#		self.read_tuple_id_width,
		#		self.genome_id_width,
		#		self.chr_id_width,
		#		self.coor_width,
		#	)

	def combine(*rnf_profiles):
		"""Combine more profiles and set their maximal values.

		Args:
			*rnf_profiles (rnftools.rnfformat.RnfProfile): RNF profile.
		"""

		for rnf_profile in rnf_profiles:
			self.prefix_width=max(self.prefix_width,rnf_profile.prefix_width)
			self.read_tuple_id_width=max(self.read_tuple_id_width,rnf_profile.read_tuple_id_width)
			self.genome_id_width=max(self.genome_id_width,rnf_profile.genome_id_width)
			self.chr_id_width=max(self.chr_id_width,rnf_profile.chr_id_width)
			self.coor_width=max(self.coor_width,rnf_profile.coor_width)

	def load(self,read_tuple_name):
		"""Load RNF values from a read tuple name.

		Args:
			read_tuple_name (str): Read tuple name which the values are taken from.
		"""
		self.prefix_width=0
		self.read_tuple_id_width=0
		self.genome_id_width=0
		self.chr_id_width=0
		self.coor_width=0

		parts=read_tuple_name.split("__")
		self.prefix_width=len(parts[0])
		self.read_tuple_id_width=len(parts[1])

		segments=parts[2][1:-1].split("),(")
		for segment in segments:
			int_widths=list(map(len,segment.split(",")))
			self.genome_id_width=max(self.genome_id_width,int_widths[0])
			self.chr_id_width=max(self.chr_id_width,int_widths[1])
			self.coor_width=max(self.coor_width,int_widths[2],int_widths[3])

	def apply(self,
				read_tuple_name,
				read_tuple_id=None,
				synchronize_widths=True
			):
		"""Apply profile on a read tuple name and update read tuple ID.

		Args:
			read_tuple_name (str): Read tuple name to be updated.
			read_tuple_id (id): New read tuple ID.
			synchronize_widths (bool): Update widths (in accordance to this profile).
		"""
		parts=read_tuple_name.split("__")
		parts[0]=self._fill_right(parts[0],"-",self.prefix_width)
		if read_tuple_id is not None:
			parts[1]="{:x}".format(read_tuple_id)
		parts[1]=self._fill_left(parts[1],"0",self.read_tuple_id_width)

		if synchronize_widths:
			new_segments=[]
			segments=parts[2][1:-1].split("),(")
			for segment in segments:
				values=segment.split(",")
				values[0]=values[0].zfill(self.genome_id_width)
				values[1]=values[1].zfill(self.chr_id_width)
				values[3]=values[3].zfill(self.coor_width)
				values[4]=values[4].zfill(self.coor_width)
				new_segments.append("("+",".join(values)+")")
			parts[2]=",".join(new_segments)

		return "__".join(parts)

	def check(self,read_tuple_name):
		"""Check if the given read tuple name satisfies this profile.

		Args:
			read_tuple_name (str): Read tuple name.
		"""

		parts=read_tuple_name.split("__")

		if len(parts[0])!= self.prefix_width or len(parts[1])!=self.read_tuple_id_width:
			return False

		segments=parts[2][1:-1].split("),(")
		for segment in segments:
			int_widths=list(map(len,segment.split(",")))
			if self.genome_id_width !=int_widths[0]:
				return False
			if self.chr_id_width != int_widths[1]:
				return False
			if self.coor_width!=int_widths[3] or self.coor_width!=int_widths[4]:
				return False

		return True

	def get_rnf_name(self,read_tuple):
		"""Get well-formatted RNF representation of a read tuple.

		read_tuple (rnftools.rnfformat.ReadTuple): Read tuple.
		"""
		return read_tuple.stringize(
					self
				)


	@staticmethod
	def _fill_left(string,character,length):
		return (length-len(string))*character + string

	@staticmethod
	def _fill_right(string,character,length):
		return string + (length-len(string))*character

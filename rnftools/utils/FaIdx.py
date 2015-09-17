class FaIdx:
	"""Class for loading FASTA indexes.

	Args:
		fai_fo (file object): FASTA index (FAI) file to be loaded.

	Attributes:
		dict_chr_ids (dict):  FASTA IDs of chromosomes (chr -> id).
		dict_ids_chr (dict):  FASTA IDs of chromosomes (id -> chr).
		dict_chr_lengths (dict): Lengths of chromosomes (chr -> length).
		number_of_chromosomes (int): Number of chromosomes in the corresponding FASTA file.
		chr_id_width (int): Length of strings representing chromosome number.
		coor_width (int): Length of string representing coordinates.
	"""

	def __init__(self, fai_fo):
		self._dict_chr_ids = {}
		self._dict_ids_chr = {}
		self._dict_chr_lengths = {}

		if fai_fo is not None:
			# parsing FAI file
			"""
			   FAI format

			1) the name of the sequence
			2) the length of the sequence
			3) the offset of the first base in the file
			4) the number of bases in each fasta line
			5) the number of bytes in each fasta line
			"""

			pairs = []

			for (i,line) in enumerate(fai_fo,start=1):
				if line.strip()!="":
					parts=line.split("\t")
					chromosome=parts[0]
					length=int(parts[1])
					pairs.append( (chromosome,length) )
			self.load_from_list(pairs)

	# pairs:list of (chromosome, length)
	def load_from_list(self, pairs):
		assert self._dict_chr_ids=={}

		for i,(chromosome, length) in enumerate(pairs,start=1):
			self._dict_chr_ids[chromosome]=i
			self._dict_ids_chr[i]=chromosome
			self._dict_chr_lengths[chromosome]=length

		self._number_of_chromosomes=len(self._dict_chr_ids)
		self._chr_id_width=len(str(self._number_of_chromosomes))
		self._coor_width=len(str(max(self._dict_chr_lengths.values())))

	@property
	def dict_chr_ids(self):
		return self._dict_chr_ids
	
	@property
	def dict_ids_chr(self):
		return self._dict_ids_chr
	
	@property
	def dict_chr_lengths(self):
		return self._dict_chr_lengths
	
	@property
	def number_of_chromosomes(self):
		return self._number_of_chromosomes
	
	@property
	def chr_id_width(self):
		return self._chr_id_width
	
	@property
	def coor_width(self):
		return self._coor_width
	
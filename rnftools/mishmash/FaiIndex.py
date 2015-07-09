

class FaiIndex:
	"""Class for loading FASTA indexes.

	Args:
		fai_file (file): FASTA index (FAI) file to be loaded.

	Attributes:
		self.dict_chr_ids (dict):  FASTA IDs of chromosomes (chr -> id).
		self.dict_chr_lengths (dict): Lengths of chromosomes (chr -> length).
		number_of_chromosomes (int): Number of chromosomes in the corresponding FASTA file.
		chr_id_width (int): Length of strings representing chromosome number.
		coor_width (int): Length of string representing coordinates.
	"""

	def __init__(self, fai_fo):
		self.dict_chr_ids = {}
		self.dict_chr_lengths = {}

		# parsing FAI file
		"""
		   FAI format

		1) the name of the sequence
		2) the length of the sequence
		3) the offset of the first base in the file
		4) the number of bases in each fasta line
		5) the number of bytes in each fasta line
		"""

		i=1
		for line in fai_fo:
			if line.strip()!="":
				parts=line.split("\t")
				chr=parts[0]
				chr_len=int(parts[1])
				self.dict_chr_ids[chr]=i
				self.dict_chr_lengths[chr]=chr_len
				i+=1

		self.number_of_chromosomes=len(self.dict_chr_ids)
		self.chr_id_width=len(str(self.number_of_chromosomes))
		self.coor_width=len(str(max(self.dict_chr_lengths.values())))



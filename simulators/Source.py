import abc

"""
	Abstract class for a source of reads
"""
class Source(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, name, fasta):
		self.name=name
		self.source_id=len(_MISHMASH_SOURCES_)+1
		_MISHMASH_SOURCES_[self.source_id] = self
		self.fasta=fasta
		self.dict_chr_ids = {}
		self.dict_chr_lengths = {}

	############################################################################
	############################################################################


	def get_dir(self):
		return name

	"""
		Get source ID compatible to the specification (without padding),
		it is assigned automatically
	"""
	def get_source_id(self):
		return self.source_id

	"""
		Get input Fasta file (registered when object was created),
		it can be an empty list
	"""
	def get_input_fa(self):
		return self.fasta

	"""
		Get required programs (with respect to snakemake-lib)
	"""
	@abc.abstractmethod
	def get_input_progs(self):
		return

	############################################################################
	############################################################################

	"""
		Get output FQ (filename)
	"""
	def get_output_fq(self):
		return single_fq_fn(self.name,self.source_id)

	"""
		Get other output files (will be marked as temp)
	"""
	@abc.abstractmethod
	def get_output_tmp(self):
		return

	############################################################################
	############################################################################

	"""
		Run read simulation
	"""
	@abc.abstractmethod
	def run(self):
		return

	def recode_sam_reads(self,sam,faidx,source=0,number_of_reads=10**9):
		self.load_fai(faidx)
		last_read_name=[]
		read_id=1

		id_str_size=len(format(number_of_reads,'x'))

		rn_formatter = RnFormatter(
				id_str_size=id_str_size,
				source_str_size=2,
				chr_str_size=self.chr_str_size,
				pos_str_size=self.pos_str_size
			)

		blocks_buffer=[]
		sequences_buffer=[]
		old_read_name=""


		cigar_reg_sh=re.compile("([0-9]+)([MDNP=X])")
		with open(sam, "r") as f:
			with open(self.get_output_fq(), "w") as ff:
				for line in f:
					line=line.strip()
					if line!="" and line[0]!="@":
						parts=line.split("\t")
						read_name=parts[0]

						if read_name!=last_read_name and read_name!="":
							read = Read(blocks=blocks_buffer)
							new_read_name = rn_formatter.process_read(read_id,read)
							read_id+=1

							for i in range(1,len(sequences_buffer)+1):
								# fixme: when only one block
								ff.write("".
									join([
										"@",new_read_name,"/",str(i) if len(sequences_buffer)>1 else "",os.linesep,
										sequences_buffer[i-1][0],os.linesep,
										"+",os.linesep,
										sequences_buffer[i-1][1],os.linesep
									]))
							last_read_name = read_name
							blocks_buffer = []
							sequences_buffer = []


						flags=int(parts[1])
						direction="F"

						#unmapped? ...skip
						if flags & 4:
							continue

						if flags & 16:
							direction="R"

						chr_id=self.dict_chr_ids[ parts[2] ] if self.dict_chr_ids!={} else "0"
						cigar=parts[5].strip()
						bases=parts[9]
						qualities=parts[10]

						left=int(parts[3])
						right=left-1
						for (steps,operation) in cigar_reg_sh.findall(cigar):
							right+=int(steps)

						block=Block(
								source=source,
								chr=chr_id,
								direction=direction,
								left=left,
								right=right
							)

						blocks_buffer.append(block)
						#fixme -- reverse complement of bases + qualities
						sequences_buffer.append( (bases,qualities) )


	"""
		name2number

	"""
	def load_fai(self,faidx):
		self.dict_chr_ids = {}
		self.dict_chr_lengths = {}

		# parsing FAI file
		with open(faidx) as f:
			"""
			   FAI format

			1) the name of the sequence
			2) the length of the sequence
			3) the offset of the first base in the file
			4) the number of bases in each fasta line
			5) the number of bytes in each fasta line
			"""

			i=1
			for line in f:
				if line.strip()!="":
					parts=line.split("\t")
					chr=parts[0]
					chr_len=int(parts[1])
					self.dict_chr_ids[chr]=i
					self.dict_chr_lengths[chr]=chr_len
					i+=1

		self.number_of_chromosomes=len(self.dict_chr_ids)
		self.chr_str_size=len(str(self.number_of_chromosomes))
		self.pos_str_size=len(str(max(self.dict_chr_lengths.values())))
		



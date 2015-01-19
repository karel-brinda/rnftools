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
		pass

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
		dict_chr=self.load_fai(faidx)
		last_read_name=[]
		read_id=1


		number_of_chromosomes=len(dict_chr)
		chr_str_size=len(str(number_of_chromosomes))
		pos_str_size=max( [ len(str(x)) for x in  dict_chr] )
		id_str_size=len(format(number_of_reads,'x'))

		rn_formatter = RnFormatter(
				id_str_size=id_str_size,
				source_str_size=2,
				chr_str_size=chr_str_size,
				pos_str_size=pos_str_size
			)

		blocks_buffer=[]
		sequences_buffer=[]
		old_read_name=""

		cigar_reg=re.compile(r"(\d+[MIDNSHP=X])+")
		cigar_reg_sh=re.compile(r"(\d+)([MDNP=X])")

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
								ff.write("".join(["@",new_read_name,"/",i,os.linesep]))
								ff.write(sequences_buffer[i-1][0],os.linesep)
								ff.write("+",os.linesep)
								ff.write(sequences_buffer[i-1][1],os.linesep)
							last_read_name = read_name
							blocks_buffer = []
							sequences_buffer = []


						flags=int(parts[1])

						#unmapped? ...skip
						if flags & 4:
							continue

						chr_number=int( dict_chr[ parts[2] ] )
						cigar=parts[5]
						bases=parts[9]
						qualities=parts[10]

						left=int(parts[3])
						right=left-1
						m=cigar_reg.search(cigar)
						for g in res.groups():
							gg=cigar_reg_sh.match(str(g))
							if gg:
								right+=int(gg.group(1))

						block=Block(
								source=source,
								chr=chr_number,
								direction=direction,
								left=left,
								right=right
							)

						blocks_buffer.append(block)
						#fixme -- reverse complement of bases + qualities
						sequences_buffer.append( (bases,qualities) )



	def load_fai(self,faidx):
		dict_chr = {}

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
					seq_len=parts[1]
					dict_chr[chr]=i
					i+=1
		return dict_chr
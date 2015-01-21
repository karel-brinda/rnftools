import os

#
# AUXILIARY FUNCTIONS
#

# TODO
# - check parameters
#
#

dwgsim_pattern = re.compile('@(.*)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_(([0-9abcdef])+)')

class DwgSim(Source):
	#TODO:estimate_unknown_values=False,
	def __init__(self,
			name,
			fasta,
			coverage=1,
			number_of_reads=0,
			read_length_1=100,
			read_length_2=0,
			other_params=""
		):
		
		super().__init__(
				name=name,
				fasta=fasta
			)

		self.read_length_1=read_length_1
		self.read_length_2=read_length_2
		#TODO: self.estimate_unknown_values=estimate_unknown_values
		self.other_params=other_params

		if coverage*number_of_reads!=0:
			raise ValueError("coverage or number_of_reads must be equal to zero")
		elif number_of_reads == 0:
			genome_size=os.stat(self.fasta).st_size
			self.number_of_reads=int(coverage*genome_size/(self.read_length_1+read_length_2))
		else:
			self.number_of_reads=number_of_reads

		self.dwg_prefix=os.path.join(
			self.get_dir(),
			self.name+"tmp.{}".format(self.source_id)
		)


	def get_input_progs(self):
		return [PROG_DWGSIM]

	def get_output_tmp(self):
		return 	[
					self.dwg_prefix+".bwa.read1.fastq",
					self.dwg_prefix+".bwa.read2.fastq",
					self.dwg_prefix+".bfast.fastq",
					self.dwg_prefix+".mutations.vcf",
					self.dwg_prefix+".mutations.txt"
			]

	def run(self):
		shell("""
				{dwgsim} -1 {rlen1} -2 {rlen2} -y 0 -N {nb} {other_params} {fasta} {pref}
			""".format(
				dwgsim=self.get_input_progs()[0],
				fasta=self.fasta,
				pref=self.dwg_prefix,
				nb=self.number_of_reads,
				rlen1=self.read_length_1,
				rlen2=self.read_length_2,
				other_params=self.other_params
			)
		)
		self.recode_dwgsim_reads(
			fastq1=self.dwg_prefix+".bfast.fastq",
			fastq2=self.get_output_fq(),
			faidx=self.fasta+".fai",
			source=self.source_id
		)
		for x in self.get_output_tmp():
			os.remove(x)


	def recode_dwgsim_reads(self,fastq1,fastq2,faidx,source=0,number_of_reads=10**9):
		max_seq_len=0

		self.load_fai(faidx)
		id_str_size=len(format(number_of_reads,'x'))

		rn_formatter = RnFormatter(
				id_str_size=id_str_size,
				source_str_size=2,
				chr_str_size=self.chr_str_size,
				pos_str_size=self.pos_str_size
			)

		#one or two ends?
		single_end=os.stat(self.dwg_prefix+".bwa.read2.fastq").st_size==0

		# parsing FQ file
		last_new_read_name=""
		with open(fastq1,"r+") as f1:
			with open(fastq2,"w+") as f2:
				i=0
				for line in f1:
					if i%4==0:
						(line1,line2,line3,line4)=("","","","")

						m = dwgsim_pattern.search(line)

						"""
						    DWGSIM read name format

						1)  contig name (chromsome name)
						2)  start end 1 (zero-based)
						3)  start end 2 (zero-based)
						4)  strand end 1 (0 - forward, 1 - reverse)
						5)  strand end 2 (0 - forward, 1 - reverse)
						6)  random read end 1 (0 - from the mutated reference, 1 - random)
						7)  random read end 2 (0 - from the mutated reference, 1 - random)
						8)  number of sequencing errors end 1 (color errors for colorspace)
						9)  number of SNPs end 1
						10) number of indels end 1
						11) number of sequencing errors end 2 (color errors for colorspace)
						12) number of SNPs end 2
						13) number of indels end 2
						14) read number (unique within a given contig/chromsome)
						"""

						contig_name     = m.group(1)
						start_1         = int(m.group(2))
						start_2         = int(m.group(3))
						direction_1     = "F" if int(m.group(4))==0 else "R"
						direction_2     = "F" if int(m.group(5))==0 else "R"
						random_1        = bool(m.group(6))
						random_2        = bool(m.group(7))
						seq_err_1       = int(m.group(8))
						snp_1           = int(m.group(9))
						indels_1        = int(m.group(10))
						seq_err_2       = int(m.group(11))
						snp_2           = int(m.group(12))
						indels_2        = int(m.group(13))
						read_id         = int(m.group(14),16)

						chr_id = self.dict_chr_ids[contig_name] if self.dict_chr_ids!={} else "0"

					elif i%4==1:
						line2=line
						read_length=len(line2.strip())

						block1=Block(
								source=source,
								chr=chr_id,
								direction=direction_1,
								left=start_1 + 1,
								right=0
							)

						block2=Block(
								source=source,
								chr=chr_id,
								direction=direction_2,
								left=start_2 + 1,
								right=0
							)
						if single_end:
							read = Read(blocks=[block1])
						else:
							read = Read(blocks=[block1,block2])
						new_read_name = rn_formatter.process_read(read_id,read)

						if single_end:
							ends_suffix=""
						elif last_new_read_name!=new_read_name:
							ends_suffix="/1"
						else:
							ends_suffix="/2"

						line1="".join(["@",new_read_name,ends_suffix,os.linesep])
						last_new_read_name=new_read_name

					elif i%4==2:
						line3=line

					else:
						line4=line

						f2.write(line1)
						f2.write(line2)
						f2.write(line3)
						f2.write(line4)

					i+=1

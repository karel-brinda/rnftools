import os

#
# AUXILIARY FUNCTIONS
#

# TODO
# - check parameters
#
#

dwgsim_pattern = re.compile('@(.*)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_(([0-9]|a|b|c|d|e|f)+)')

class DwgSim(Source):
	def __init__(self,
			name,
			fasta,
			coverage=1,
			read_count=0,
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
		self.other_params=other_params

		self.dwg_prefix=os.path.join(
			self.get_dir(),
			self.name+"___tmp___.{}".format(self.source_id)
		)


		if coverage*read_count!=0:
			# todo: error - only 1 can be set
			pass
		else:
			if read_count != 0:
				self.read_count=read_count
			else:
				genome_size=os.stat(self.fasta).st_size
				print(os.stat(self.fasta).st_size)
				self.read_count=int(coverage*genome_size/(self.read_length_1+read_length_2))

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
				nb=self.read_count,
				rlen1=self.read_length_1,
				rlen2=self.read_length_2,
				other_params=self.other_params
			)
		)
		self.recode_dwgsim_reads(
			fastq1=self.dwg_prefix+".bwa.read1.fastq",
			fastq2=self.get_output_fq(),
			faidx=self.fasta+".fai",
			source=self.source_id
		)
		for x in self.get_output_tmp():
			os.remove(x)


	def recode_dwgsim_reads(self,fastq1,fastq2,faidx,source=0,number_of_reads=10**9):
		max_seq_len=0

		dict_chr=self.load_fai(faidx)
		number_of_chromosomes=len(dict_chr)
		chr_str_size=len(str(number_of_chromosomes))
		pos_str_size=max( [ len(str(x)) for x in  dict_chr] )
		id_str_size=len(format(number_of_reads,'x'))

		# parsing FQ file
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

						contig_name 	=	 m.group(1)
						start_end_1 	=	 int(m.group(2))
						start_end_2 	=	 int(m.group(3))
						strand_end_1 	=	 "F" if int(m.group(4))==0 else "R"
						strand_end_2 	=	 "F" if int(m.group(5))==0 else "R"
						random_1 		=	 bool(m.group(6))
						random_2 		=	 bool(m.group(7))
						seq_err_1 		=	 int(m.group(8))
						snp_1		 	=	 int(m.group(9))
						indels_ 		=	 int(m.group(10))
						seq_err_2 		=	 int(m.group(11))
						snp_2 			=	 int(m.group(12))
						indels_2 		=	 int(m.group(13))
						read_number 	=	 int(m.group(14),16)

						chr = dict_chr[contig_name] if dict_chr!={} else "0"

					elif i%4==1:
						line2=line
						read_length=len(line2.strip())
						starting_pos=start_end_1 + 1
						end_pos=starting_pos+read_length-1
						block="({},{},{},{},{})".format(
								source,
								str(chr).zfill(chr_str_size),
								strand_end_1,
								str(starting_pos).zfill(pos_str_size),
								str(end_pos).zfill(pos_str_size)
							)
						blocks=[block]

						line1="@{}__{}__{}__{}".format(
								"",
								format(read_number,'x').zfill(id_str_size),
								"".join(blocks),
								os.linesep
							)

					elif i%4==2:
						line3=line

					else:
						line4=line

						f2.write(line1)
						f2.write(line2)
						f2.write(line3)
						f2.write(line4)

					i+=1

#
# AUXILIARY FUNCTIONS
#

dwgsim_pattern = re.compile('@(.*)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_(([0-9]|a|b|c|d|e|f)+)')

def dwgsim_output_files(name,batch_id):
	prefix = "tmp/"+single_fq_fn(name,batch_id,"dwgsim")
	return 	[
				prefix+".bfast.fastq",
				prefix+".bwa.read1.fastq",
				prefix+".bwa.read2.fastq",
				prefix+".mutations.vcf",
				prefix+".mutations.txt"
			]		

def recode_dwgsim_reads(fastq1,fastq2,faidx="",source=0,number_of_reads=10**9):
	dict_chr = {}
	max_seq_len=0

	# parsing FAI file
	if faidx!="":
		with open(faidx) as f:
			## FAI format
			# the name of the sequence
			# the length of the sequence
			# the offset of the first base in the file
			# the number of bases in each fasta line
			# the number of bytes in each fasta line

			i=1
			for line in f:
				if line.strip()!="":
					parts=line.split("\t")
					chr=parts[0]
					seq_len=parts[1]
					max_seq_len=max(max_seq_len,int(seq_len))
					dict_chr[chr]=i
					i+=1
	#print(dict_chr)
	number_of_chromosomes=len(dict_chr)
	chr_str_size=len(str(number_of_chromosomes))
	pos_str_size=len(str(max_seq_len))
	id_str_size=len(format(number_of_reads,'x'))

	# parsing FQ file
	with open(fastq1,"r+") as f1:
		with open(fastq2,"w+") as f2:
			i=0
			for line in f1:
				if i%4==0:
					(line1,line2,line3,line4)=("","","","")

					m = dwgsim_pattern.search(line)

					## DWGSIM read name format
					#contig name (chromsome name)
					#start end 1 (zero-based)
					#start end 2 (zero-based)
					#strand end 1 (0 - forward, 1 - reverse)
					#strand end 2 (0 - forward, 1 - reverse)
					#random read end 1 (0 - from the mutated reference, 1 - random)
					#random read end 2 (0 - from the mutated reference, 1 - random)
					#number of sequencing errors end 1 (color errors for colorspace)
					#number of SNPs end 1
					#number of indels end 1
					#number of sequencing errors end 2 (color errors for colorspace)
					#number of SNPs end 2
					#number of indels end 2
					#read number (unique within a given contig/chromsome)

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
import time
import re

include: "snakemake-lib/progs.py"

# misto slovniku objekt
# chromosome name => number


def single_fq_fn(name,batch_id,simulator):
	return "{}.part_{}.{}.fq".format(name,batch_id,simulator)

def joined_fq_fn(name):
	return "{}.joined.fq".format(name)

# todo: check simulators

rule joined_fq:
	input:
		[single_fq_fn("{name}",i,samples[i]["simulator"]) for i in range(len(samples))]
	output:
		joined_fq_fn("{name}")
	run:
		shell("cat /dev/null > {output[0]}")
		for i in range(len(input)):
			shell("cat " + input[0] + " >> {output[0]}")

def prog_from_progname(wildcards):
	print(wildcards)
	if samples[int(wildcards.batch_id)]["simulator"]=="DWGSIM":
		return DWGSIM

def fasta_from_id(wildcards):
	return samples[int(wildcards.batch_id)]["fasta"]

def dwgsim_output_files(name,batch_id):
	prefix = "tmp/"+single_fq_fn(name,batch_id,"dwgsim")
	return 	[
				prefix+".bfast.fastq",
				prefix+".bwa.read1.fastq",
				prefix+".bwa.read2.fastq",
				prefix+".mutations.vcf",
				prefix+".mutations.txt"
			]		

def recode_dwgsim_reads(fastq1,fastq2,source=0):
	with open(fastq1,"r+") as f1:
		with open(fastq2,"w+") as f2:
			i=0
			for line in f1:
				if i%4==0:
					#@gi|561108321|ref|NC_018143.2|_2359637_1_0_1_0_0_1:0:0_0:0:0_23e70
					reg='@(.*)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_([0-9]+):([0-9]+):([0-9]+)_(([0-9]|a|b|c|d|e|f)+)'
					m = re.search(reg, line)

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

					block="({},{},{},{},{})".format(source,contig_name,strand_end_1,start_end_1,0)

					print (block, read_number)
					#print("\t".join( m.groups()[1:] ))
					#print(line)

				else:
					f2.write(line)
				i+=1


rule sub_fq_dwgsim:
	output:
		temp(single_fq_fn("{name}","{batch_id}","dwgsim")),
		#single_fq_fn("{name}","{batch_id}","dwgsim")+".bfast.fastq",
		temp(dwgsim_output_files("{name}","{batch_id}"))
	input:
		DWGSIM,
		fasta_from_id
	params:
		DWGSIM=DWGSIM
	run:
		shell(
			"""
				{params.DWGSIM} -1 30 -2 0 -y 0 -C 1 {input[1]} tmp/{output[0]}
			"""
		)
		recode_dwgsim_reads(output[1],output[0],source=int(wildcards.batch_id))

#rule:
#    input: 
#        ART_ILLUMINA,
#        EXAMPLE_FASTA
#    shell:
#        "art_illumina -sam -i {input[1]} -p -l 50 -f 20 -m 200 -s 10 -o paired_dat"

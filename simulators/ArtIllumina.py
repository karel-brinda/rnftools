import os
import pysam

#
# AUXILIARY FUNCTIONS
#

# TODO
# - check parameters
#
#

class ArtIllumina(Source):
	def __init__(self,
			name,
			fasta,
			coverage=1,
			number_of_reads=0,
			read_length=100,
			other_params=""
		):
		
		super().__init__(
				name=name,
				fasta=fasta
			)

		#fixme pairs
		if coverage*number_of_reads!=0:
			raise ValueError("coverage or number_of_reads must be equal to zero")
		elif number_of_reads == 0:
			self.coverage=coverage
		else:
			genome_size=os.stat(self.fasta).st_size
			self.coverage=1.0*number_of_reads*read_length/genome_size


		self.read_length=read_length
		self.other_params=other_params

		self.art_prefix=os.path.join(
			self.get_dir(),
			self.name+"tmp.{}".format(self.source_id)
		)

	def get_input_progs(self):
		return [PROG_ART_ILLUMINA]

	def get_output_tmp(self):
		return 	[
					self.art_prefix+".bwa.read1.fastq",
					self.art_prefix+".bwa.read2.fastq",
					self.art_prefix+".bfast.fastq",
					self.art_prefix+".mutations.vcf",
					self.art_prefix+".mutations.txt"
			]


	def run(self):
		shell("""
				{art_il} -sam -na \
					-i {fasta} \
					-l {rlen} \
					{other_params} \
					-f {coverage} \
					-p -m 200 -s 10 \
					 -o {o_pref}
			""".format(
				art_il=self.get_input_progs()[0],
				fasta=self.fasta,
				rlen=self.read_length,
				other_params=self.other_params,
				coverage=self.coverage,
				o_pref=self.art_prefix
			)
		)
		self.recode_sam_reads(
			sam=self.art_prefix+".sam",
			faidx=self.fasta+".fai",
			source=self.source_id
		)
#		for x in self.get_output_tmp():
#			os.remove(x)



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
			read_count=0,
			read_length=100,
			other_params=""
		):
		
		super().__init__(
				name=name,
				fasta=fasta
			)
		self.read_length=read_length
		self.other_params=other_params

		self.art_prefix=os.path.join(
			self.get_dir(),
			self.name+"___tmp___.{}".format(self.source_id)
		)

		"art_illumina -sam -i reference.fa -l 50 -f 10 -o single_dat"
		"art_illumina -i example_fasta.fa -l 100  -sam -na -f 10 -o testtt"

		if coverage*read_count!=0:
			# todo: error - only 1 can be set
			pass
		else:
			if read_count != 0:
				self.read_count=read_count
			else:
				genome_size=os.stat(self.fasta).st_size
				print(os.stat(self.fasta).st_size)
				#FIXME: pairs
				self.read_count=int(coverage*genome_size/(self.read_length))

	def get_input_progs(self):
		return [PROG_ART_ILLUMINA]

	def get_output_tmp(self):
		return 	[					
			]

	def run(self):
		shell("""
				{art_il} -sam -na \
					-i {fasta} \
					-l {rlen} \
					{other_params} \
					 -f 0.1 -p -m 200 -s 10 \
					 -o {o_pref}
			""".format(
				art_il=self.get_input_progs()[0],
				fasta=self.fasta,
				rlen=self.read_length,
				other_params=self.other_params,
				o_pref=self.art_prefix			)
		)
		self.recode_sam_reads(
			sam=self.art_prefix+".sam",
			faidx=self.fasta+".fai",
			source=self.source_id
		)
#		for x in self.get_output_tmp():
#			os.remove(x)



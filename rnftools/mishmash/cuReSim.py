from rnftools import mishmash

import os
import smbl
import snakemake
import re

#
# AUXILIARY FUNCTIONS
#

class CuReSim(mishmash.Source):
	def __init__(self,
			fa,
			coverage=0,
			number_of_reads=0,
			read_length_1=100,
			read_length_2=0,
			other_params="",
			rng_seed=1
		):
		
		if read_length_2!=0:
			raise ValueError("CuReSim supports only single-end reads")

		super().__init__(
				fa=fa,
				ends=1,
				rng_seed=rng_seed,
				number_of_threads=9999,
			)


		self.read_length_1=read_length_1
		self.read_length_2=read_length_2
		self.other_params=other_params

		self.number_of_reads=number_of_reads
		self.coverage=coverage


	def get_input(self):
		return [
				smbl.prog.CURESIM,
				self._fa_fn,
				self._fai_fn,
			]

	def get_output(self):
		return [
				self._fq_fn,
#				os.path.join(
#					self.get_dir(),
#					"output.fastq",
#				),
#				os.path.join(
#					self.get_dir(),
#					"log.txt",
#				),
			]

	# TODO: find out how it is with RNG seeds
	def create_fq(self):

		if self.number_of_reads == 0:
			genome_size=os.stat(self._fa_fn).st_size
			self.number_of_reads=int(self.coverage*genome_size/(self.read_length_1+self.read_length_2))

		snakemake.shell("""
				cd {dir}
				java -Xmx8g -jar \
				{curesim} \
				-f "{fa}" \
				-n {nb} \
				-m {rlen1} \
				-r 0 \
				-sd 0 \
				-y 0 \
				{other_params} \
				> /dev/null
			""".format(
				dir=self.get_dir(),
				curesim=smbl.prog.CURESIM,
				fa=self._fa_fn,
				nb=self.number_of_reads,
				rlen1=self.read_length_1,
				other_params=self.other_params,
				rng_seed=self._rng_seed,
			)
		)
		self.recode_curesim_reads(
			os.path.join(
					self.get_dir(),
					"output.fastq",
				)
		)

	def recode_curesim_reads(self,old_fq,number_of_reads=10**9):
		curesim_pattern = re.compile('@(.*)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)')
		"""
			CuReSim read name format

			@<#1>_<#2>_<#3>_<#4>_<#5>_<#6>_<#7>_<#8>

			1: contig name
			2: original position
			3: strand (0=forward;1=reverse)
			4: random read (0=non-random;1=random)
			5: number of insertions
			6: number of deletions
			7: number of substitution
			8: read number (unique within a genome)
		"""

		max_seq_len=0

		self.load_fai()
		id_str_size=len(format(number_of_reads,'x'))

		rn_formatter = smbl.RnFormatter(
				id_str_size=id_str_size,
				source_str_size=2,
				chr_str_size=self.chr_str_size,
				pos_str_size=self.pos_str_size
			)

		# parsing FQ file
		last_new_read_name=""
		read_id=0
		with open(old_fq,"r+") as f1:
			with open(self._fq_fn,"w+") as f2:
				i=0
				for line in f1:
					if i%4==0:
						(line1,line2,line3,line4)=("","","","")

						m = curesim_pattern.search(line)
						if m is None:
							raise ValueError("Read '{}' was not generated by CuReSim.".format(line[1:]))

						contig_name     = m.group(1)
						start_pos       = int(m.group(2))
						direction       = "R" if int(m.group(3)) else "F"
						random          = bool(m.group(4))
						ins_nb          = int(m.group(5))
						del_nb          = int(m.group(6))
						subst_nb        = int(m.group(7))
						rd_id           = int(m.group(8))

						end_pos         = start_pos - 1 - ins_nb + del_nb

						chr_id=0
						#TODO: uncomment when the chromosome naming bug in curesim is corrected
						#chr_id = self.dict_chr_ids[contig_name] if self.dict_chr_ids!={} else "0"

					elif i%4==1:
						line2=line
						read_length=len(line2.strip())
						end_pos += read_length

						block=smbl.Block(
								source=self.source_id,
								chr=chr_id,
								direction=direction,
								left=start_pos + 1,
								right=end_pos,
							)

						read = smbl.Read(blocks=[block],read_id=read_id+1,suffix="[single-end,curesim]")
						new_read_name = rn_formatter.process_read(read)

						read_id+=1

						line1="".join(["@",new_read_name,os.linesep])
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

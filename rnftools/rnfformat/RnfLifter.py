from .chain_dict import ChainDict
import rnftools

import re

class RnfLifter:

	def __init__(self,
				chain_fn,
				fai_fn=,
			):		if 
		self._fai_fn=fai_fn
		self._chain_fn=chain_fn

		self._chain_dict=ChainDict(chain_fn)
		self._fai_index=rnftools.mishmash.FaiIndex(fai_fn)

		self._reg_block=re.compile(r"(\(([0-9]+),[0-9]+,[FRN],([0-9]+),([0-9]+)\))")

	def lift_rnf_name(self,rnf_name):
		for occur in self._reg_block.finditer(rnf_name):
			print(occur)
		return rnf_name

	def lift_fastq(self, fastq_in_fo, fastq_out_fo):
		for i, line in enumerate(fastq_in_fo,start=0):
			if i%4==0:
				fastq_out_fo.write(self.lift_rnf_name(line))
			else
				fastq_out_fo.write(line)
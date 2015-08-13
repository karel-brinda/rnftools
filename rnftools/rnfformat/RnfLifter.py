import rnftools
import pysam

import re

class RnfLifter:

	def __init__(self,
				chain_fn,
				fai_fn,
			):
		self._fai_fn=fai_fn
		self._chain_fn=chain_fn

		with open(chain_fn) as chain_fo:
			self._chain=rnftools.utils.Chain(chain_fo=chain_fo)
		with open(fai_fn) as fai_fo:
			self._fai_index=rnftools.utils.FaiIndex(fai_fo=fai_fo)

		self._reg_block=re.compile(r"(\(([0-9]+),[0-9]+,[FRN],([0-9]+),([0-9]+)\))")

	def lift_rnf_name(self,rnf_name):
		for occur in self._reg_block.finditer(rnf_name):
			groups=occur.groups()
			chrom_id=int(groups[1])
			chrom=self._fai_index.dict_ids_chr[chrom_id]
			o_left=groups[2]
			o_right=groups[3]
			left=int(o_left)
			right=int(o_right)
			if left!=0:
				n_left=self._chain.one_based_transl(chrom,left)
				f_new_left=str(n_left).zfill(len(o_left))
				rnf_name=rnf_name.replace(",{},".format(o_left),",{},".format(f_new_left))
			if right!=0:
				new_right=self._chain.one_based_transl(chrom,right)
				f_new_right=str(n_right).zfill(len(o_right))
				rnf_name=rnf_name.replace(",{})".format(o_right),",{})".format(f_new_right))
		return rnf_name

	def lift_fastq(self,
				fastq_in_fo,
				fastq_out_fo
			):
		for i, line in enumerate(fastq_in_fo,start=0):
			if i%4==0:
				fastq_out_fo.write(self.lift_rnf_name(line))
			else:
				fastq_out_fo.write(line)

	def lift_sam(self,
				sam_in_fn,
				sam_out_fn
			):
		infile = pysam.AlignmentFile(sam_in_fn, "rb" if sam_in_fn[-4:]==".bam" else "r")
		outfile = pysam.AlignmentFile(sam_out_fn, "wb" if sam_out_fn[-4:]==".bam" else "wh", template=infile)
		for s in infile:
			s.qname=self.lift_rnf_name(s.qname)
			outfile.write(s)
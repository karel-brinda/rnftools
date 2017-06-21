import rnftools
import pysam

import re, sys


class RnfLifter:
	def __init__(
			self,
			chain_fn,
			fai_fn,
			invert=False,
	):
		self._fai_fn = fai_fn
		self._chain_fn = chain_fn
		self._invert = invert

		if chain_fn is not None:
			with open(chain_fn) as chain_fo:
				self._chain = rnftools.utils.Chain(chain_fo=chain_fo, invert=invert)
		else:
			self._chain = None

		if fai_fn is not None:
			with open(fai_fn) as fai_fo:
				self._fai_index = rnftools.utils.FaIdx(fai_fo=fai_fo)
		else:
			if self._chain is not None:
				self._fai_index = self._chain.get_fasta_index()
			else:
				self._fai_index = None

		self._reg_block = re.compile(r"(\(([0-9]+),([0-9]+),[FRN],([0-9]+),([0-9]+)\))")

	def lift_rnf_name(self, rnf_name, genome_id):
		if self._chain is None:
			return rnf_name
		for occur in self._reg_block.finditer(rnf_name):
			groups = occur.groups()
			p_genome_id = int(groups[1])
			chrom_id = int(groups[2])

			# is this segment from a genome to be transformed?1
			if int(genome_id) == int(p_genome_id):
				# print("going to transform",file=sys.stderr)
				chrom = self._fai_index.dict_ids_chr[chrom_id]
				o_left = groups[3]
				o_right = groups[4]
				left = int(o_left)
				right = int(o_right)
				if left != 0:
					n_left = self._chain.one_based_transl(chrom, left)
					f_new_left = str(n_left).zfill(len(o_left))
					rnf_name = rnf_name.replace(",{},".format(o_left), ",{},".format(f_new_left))
				# print("left",o_left,"=>",f_new_left,file=sys.stderr)
				if right != 0:
					new_right = self._chain.one_based_transl(chrom, right)
					f_new_right = str(n_right).zfill(len(o_right))
					rnf_name = rnf_name.replace(",{})".format(o_right), ",{})".format(f_new_right))
				# print(o_right,"=>",f_new_right,file=sys.stderr)
		return rnf_name

	def lift_fastq(self,
			fastq_in_fo,
			fastq_out_fo,
			genome_id,
	):
		for i, line in enumerate(fastq_in_fo, start=0):
			if i % 4 == 0:
				fastq_out_fo.write(self.lift_rnf_name(line, genome_id=int(genome_id)))
			else:
				fastq_out_fo.write(line)

	def lift_sam(
			self,
			genome_id,
			sam_in_fn=None,
			bam_in_fn=None,
			sam_out_fn=None,
			bam_out_fn=None,
	):
		assert sam_in_fn is not None or bam_in_fn is not None
		assert sam_in_fn is None or bam_in_fn is None
		assert sam_out_fn is not None or bam_out_fn is not None
		assert sam_out_fn is None or bam_out_fn is None

		if sam_in_fn is None:
			in_mode = "rb"
			in_fn = bam_in_fn
		else:
			in_mode = "r"
			in_fn = sam_in_fn

		if sam_out_fn is None:
			out_mode = "wb"
			out_fn = bam_out_fn
		else:
			out_mode = "wh"
			out_fn = sam_out_fn

		infile = pysam.AlignmentFile(in_fn, in_mode)
		outfile = pysam.AlignmentFile(out_fn, out_mode, template=infile)
		for s in infile:
			s.qname = self.lift_rnf_name(s.qname, genome_id=int(genome_id))
			outfile.write(s)

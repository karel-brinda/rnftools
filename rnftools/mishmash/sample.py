import rnftools.mishmash

import smbl
import snakemake
import abc
import re
import os

class Sample:
	def __init__(self, name, ends, paired_end_mode="bwa"):
		self._name=name
		self._sources=[]
		self._dir=name
		self._ends=ends
		self._pair_end_mode=paired_end_mode

		if self._ends==1:
			self._fq_fns=[self._name+".fq"]
			self._mode="single-end"
		elif self._ends==2 and self._pair_end_mode=="bwa":
			self._fq_fns=[self._name+".1.fq",self._name+".2.fq"]
			self._mode="pair-end-bwa"
		elif self/_ends==2 and self._pair_end_mode=="bfast":
			self._fq_fns=[self._name+".fq"]
			self._mode="pair-end-bfast"

		rnftools.mishmash.add_sample(self)

		if paired_end_mode not in ["bwa","bfast"]:
			raise ValueError("paired_end_mode must be 'bwa' or 'bfast'")

	def get_name(self):
		return self._name

	def get_dir(self):
		return self._dir

	def get_sources(self):
		return self._sources

	def get_input(self):
		return [source.fq_fn() for source in self._sources]

	def add_source(self,source):
		if self._ends!=source.get_number_of_ends():
			raise ValueError(
				"It is not possible to combine reads with different number of ends in a single sample. "
				"Details: name='{}', old ends='{}', new ends='{}', source='{}'.".format(self._name,self._ends,source.get_number_of_ends(),source)
			)
		self._sources.append(source)

	def clean(self):
		for x in self._fq_fns+[self._dir]:
			snakemake.shell('rm -fR "{}"'.format(x))

	######################################
	######################################

	def fq_fns(self):
		return self._fq_fns

	######################################
	######################################

	def create_fq(self):
		# fixme: use absolute path to this file
		snakemake.shell("""
				rnf-join-fq.py \
					-i {input_fqs} \
					-m {mode} \
					-o "{output_prefix}"\
			""".format(
					input_fqs=" ".join(['"{}"'.format(source.fq_fn()) for source in self._sources]),
					mode=self._mode,
					output_prefix=self._name,
		))

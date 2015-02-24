#! /usr/bin/env python3

import argparse
import os
import random

# number of reads taken in a single run
READS_IN_GROUP=10

allowed_modes = [
			"single-end",
			"pair-end-bwa",
			"pair-end-bfast",
		]

class Mixer:
	def __init__(self,mode,input_files,output_prefix):
		self.mode=mode
		self.input_files=input_files
		self.output_prefix=output_prefix
		self.rng=random.Random()
		self.rng.seed(1)

		if args.m=="single-end":
			self.output=Output(fn_1="{}.fq".format(output_prefix),ends=1)
			self._ends=1
		elif args.m=="pair-end-bwa":
			self.output=Output(fn_1="{}.1.fq".format(output_prefix),fn_2="{}.2.fq".format(output_prefix),ends=2)
			self._ends=2
		elif args.m=="pair-end-bfast":
			self.output=Output(fn_1="{}.fq".format(output_prefix),ends=2)
			self._ends=2
		else:
			raise ValueError("Unknown mode '{}'".format(args.m))

		self.i_files=[open(fn) for fn in input_files]
		self.i_files_proc=[os.path.getsize(fn) for fn in input_files]
		self.i_files_proc=[int((100.0*x)/sum(self.i_files_proc)) for x in self.i_files_proc]
		self.i_files_weighted=[]
		for i in range(len(self.i_files)):
			self.i_files_weighted.extend(self.i_files_proc[i]*[self.i_files[i]])

	def run(self):
		while len(self.i_files_weighted)>0:
			file_id=self.rng.randint(0,len(self.i_files_weighted)-1)
			for i in range(READS_IN_GROUP*self._ends):
				if self.i_files_weighted[file_id].closed:
					del self.i_files_weighted[file_id]
					break

				ln1=self.i_files_weighted[file_id].readline()
				ln2=self.i_files_weighted[file_id].readline()
				ln3=self.i_files_weighted[file_id].readline()
				ln4=self.i_files_weighted[file_id].readline()

				if not ln1:
					self.i_files_weighted[file_id].close()
					del self.i_files_weighted[file_id]
					break
				self.output.save_read(ln1,ln2,ln3,ln4)


class Output:
	def __init__(self,ends,fn_1,fn_2=None):
		self.ends=ends
		self.fs=[open(fn_1,"w+")]
		if fn_2 is not None:
			self.fs.append(open(fn_2,"w+"))

	def __del__(self):
		for f in self.fs:
			f.close()

	def save_read(self,ln1,ln2,ln3,ln4):
		[ln1,ln2,ln3,ln4]=[ln1.strip(),ln2.strip(),ln3.strip(),ln4.strip()]
	
		if self.ends==1:
			file_id=0
			if ln1[-2]=="/":
				raise ValueError("Wrong read name '{}'. Single end read should not contain '/'.".format(ln1[1:]))
	
		else:
			if ln1[-2]!="/":
				raise ValueError("Wrong read name '{}'. A read with two ends should contain '/'.".format(ln1[1:]))
			if len(self.fs)==1:
				ln1=ln1[:-2]
				file_id=0
			else:
				if ln1[-1]=="1":
					file_id=0
				elif ln1[-1]=="2":
					file_id=1
				else:
					raise ValueError("Wrong read name '{}'.".format(ln1[1:]))

		self.fs[file_id].write("".join([ln1,os.linesep,ln2,os.linesep,ln3,os.linesep,ln4,os.linesep]))


parser = argparse.ArgumentParser()

parser.add_argument(
		'-i',
		required=True,
		metavar='inp',
		nargs='+',
		help='input FASTQ files',
	)

parser.add_argument(
		'-m',
		required=True,
		metavar='mode',
		choices=allowed_modes,
		#type=lambda x: is_valid_mode(parser,x),
		help='mode',
	)


parser.add_argument(
		'-o',
		metavar='out',
		required=True,
		help='output prefix',
	)

args = parser.parse_args()

outpref=args.o
inp_fastqs=args.i

mixer=Mixer(
		mode=args.m,
		input_files=args.i,
		output_prefix=args.o,
	)
mixer.run()


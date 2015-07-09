import os
import sys
import random
import math

# number of reads taken in a single run
READS_IN_GROUP=10

ALLOWED_MODES = [
			"single-end",
			"paired-end-bwa",
			"paired-end-bfast",
		]

class FqMerger:
	def __init__(self,mode,input_files,output_prefix):
		self.mode=mode
		self.input_files=input_files
		self.output_prefix=output_prefix
		self.rng=random.Random()
		self.rng.seed(1)

		self.i_files=[open(fn) for fn in input_files]
		self.i_files_sizes=[os.path.getsize(fn) for fn in input_files]
		self.i_files_proc=[int((100.0*x)/sum(self.i_files_sizes)) for x in self.i_files_sizes]
		self.i_files_weighted=[]
		for i in range(len(self.i_files)):
			self.i_files_weighted.extend(self.i_files_proc[i]*[self.i_files[i]])

		read_id_length_est=math.ceil(
					math.log(
						sum(self.i_files_sizes)/20,
						16,
					)
				)

		if mode=="single-end":
			self.output_files=[
				"{}.fq".format(output_prefix)
			]
			self.output=FqMergerOutput(fn_1=self.output_files[0],reads_in_tuple=1,read_id_length=read_id_length_est)
			self._reads_in_tuple=1
		elif mode=="paired-end-bwa":
			self.output_files=[
				"{}.1.fq".format(output_prefix),
				"{}.2.fq".format(output_prefix),
			]
			self.output=FqMergerOutput(fn_1=self.output_files[0],fn_2=self.output_files[1],reads_in_tuple=2,read_id_length=read_id_length_est)
			self._reads_in_tuple=2
		elif mode=="paired-end-bfast":
			self.output_files=[
				"{}.fq".format(output_prefix)
			]
			self.output=FqMergerOutput(fn_1=self.output_files[0],reads_in_tuple=2,read_id_length=read_id_length_est)
			self._reads_in_tuple=2
		else:
			raise ValueError("Unknown mode '{}'".format(mode))

	def run(self):
		print("",file=sys.stderr)
		print("Going to merge/convert RNF-FASTQ files.",file=sys.stderr)
		print("",file=sys.stderr)
		print("   mode:          ", self.mode,file=sys.stderr)
		print("   input files:   ", ", ".join(self.input_files),file=sys.stderr)
		print("   output files:  ", ", ".join(self.output_files),file=sys.stderr)
		print("",file=sys.stderr)


		while len(self.i_files_weighted)>0:
			file_id=self.rng.randint(0,len(self.i_files_weighted)-1)
			for i in range(READS_IN_GROUP*self._reads_in_tuple):
				if self.i_files_weighted[file_id].closed:
					del self.i_files_weighted[file_id]
					break

				ln1=self.i_files_weighted[file_id].readline()
				ln2=self.i_files_weighted[file_id].readline()
				ln3=self.i_files_weighted[file_id].readline()
				ln4=self.i_files_weighted[file_id].readline()

				if not ln1 or not ln2 or not ln3 or not ln4:
					self.i_files_weighted[file_id].close()
					del self.i_files_weighted[file_id]
					break
				assert ln1[0]=="@", ln1
				assert ln3[0]=="+", ln3
				self.output.save_read(ln1,ln2,ln3,ln4)
		self.output.close()


class FqMergerOutput:
	def __init__(self,reads_in_tuple,fn_1,fn_2=None,read_id_length=6):
		self.reads_in_tuple=reads_in_tuple
		self.read_id_length=read_id_length
		self.fs=[open(fn_1,"w+")]
		if fn_2 is not None:
			self.fs.append(open(fn_2,"w+"))
		self.read_tuple_counter=0

	def close(self):
		for f in self.fs:
			f.close()

	def save_read(self,ln1,ln2,ln3,ln4):
		[ln1,ln2,ln3,ln4]=[ln1.strip(),ln2.strip(),ln3.strip(),ln4.strip()]

		ln1_parts=ln1.split("__")
		ln1_parts[1]="{:x}".format(self.read_tuple_counter).zfill(self.read_id_length)
		ln1="__".join(ln1_parts)
	
		if self.reads_in_tuple==1:
			file_id=0
			if ln1[-2]=="/":
				raise ValueError("Wrong read name '{}'. Single end read should not contain '/'.".format(ln1[1:]))
			self.read_tuple_counter+=1
	
		else:
			if ln1[-2]!="/":
				raise ValueError("Wrong read name '{}'. A read with two ends should contain '/'.".format(ln1[1:]))
			if len(self.fs)==1:
				ln1=ln1[:-2]
				file_id=0
				self.read_tuple_counter+=1
			else:
				if ln1[-1]=="1":
					file_id=0
				elif ln1[-1]=="2":
					file_id=1
					self.read_tuple_counter+=1
				else:
					raise ValueError("Wrong read name '{}'.".format(ln1[1:]))

		self.fs[file_id].write("".join([ln1,os.linesep,ln2,os.linesep,ln3,os.linesep,ln4,os.linesep]))

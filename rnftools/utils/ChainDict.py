import numpy

class ChainDict:

	def __init__(self,
				chain_fn,
			):

		self._dict_array={}
		self._dict_properties={}

		self._chain_fn=chain_fn
		chrom=None
		with open(chain_fn) as f:
			while 1==1:
				line=f.readline()
				if line=="":
					break

				head=line.strip()
				if head=="":
					continue

				#print(head)
				parts=head.split()

				[
					_,
					score,
					tName,
					tSize,
					tStrand,
					tStart,
					tEnd,
					qName,
					qSize,
					qStrand,
					qStart,
					qEnd,
					iid,
				]=parts

				score =int (score)
				tSize = int(tSize)
				tStart = int(tStart)
				tEnd = int(tEnd)
				qSize = int(qSize)
				qStart = int(qStart)
				qEnd = int(qEnd)
				iid = int(iid)

				self._dict_array[tName]=numpy.empty(tSize, dtype=numpy.int32)
				self._dict_properties[tName]={
					"score":score,
					"tName":tName,
					"tSize":tSize,
					"tStrand":tStrand,
					"tStart":tStart,
					"tEnd":tEnd,
					"qName":qName,
					"qSize":qSize,
					"qStrand":qStrand,
					"qStart":qStart,
					"qEnd":qEnd,
					"id":iid,
				}

				leftp = 0
				rightp= 0

				#line=f.readline()
				#if line is None:
				#	break
				#parts=line.strip().split()
				parts=[0, tStart, qStart]
				while len(parts)>0:
					if len(parts)==1:
						parts=[int(parts[0]), tSize-tEnd, qSize-qEnd]

					if len(parts)==3:
						[size, dt, dq] = map(int, parts)
						for i in range(size):
							self._dict_array[tName][leftp]=rightp
							leftp+=1
							rightp+=1
						if dt==0 and dq>0:

							rightp+=dq

						elif dt>0:
							slope=1.0*dq/dt
							for i in range(dt):
								self._dict_array[tName][leftp]=rightp+int(round(rightp + (i*slope)))
								leftp+=1
							rightp+=dq

					#print(leftp)
					parts=f.readline().strip().split()
				assert leftp==tSize, "{} {} vs. {} {} (diff: {} {})".format(leftp,rightp,tSize,qSize,tSize-leftp,qSize-rightp)
				assert rightp==qSize, "{} {} vs. {} {} (diff: {} {})".format(leftp,rightp,tSize,qSize,tSize-leftp,qSize-rightp)


	def one_based_transl(self, chromosome, coordinate):
		assert 0<=coordinate
		#assert coordinate<=self._dict_properties[chromosome][tSize]
		if coordinate==0:
			return 0
		else:
			return self._dict_array[chromosome][coordinate-1]+1

import collections
from .ChainSequence import ChainSequence

class Chain:

	def __init__(self,
				chain_fo,
				sampling_step=10000,
				inverted=False
			):

		self._chain_fo=chain_fo
		self._inverted=inverted
		self._sampling_step=sampling_step

		tmp_chain_sequences=[]

		for line in chain_fo:
			line=line.strip()

			if line=="":
				continue

			parts=line.split()

			assert len(parts) in [1,3,13]

			if len(parts)==13:
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

				intervals_pairs = [((0,tStart),(0,qStart))]

			elif len(parts)==3:
				[size, dt, dq] = map(int, parts)
				interval_pairs.append(
						(
							(interval_pairs[-1][0][1],interval_pairs[-1][0][1]+size),
							(interval_pairs[-1][1][1],interval_pairs[-1][1][1]+size),
						)
					)
				interval_pairs.append(
						(
							(interval_pairs[-1][0][1],interval_pairs[-1][0][1]+dt),
							(interval_pairs[-1][1][1],interval_pairs[-1][1][1]+dq),
						)
					)
			elif len(parts)==1:
				[size] = map(int, parts)
				interval_pairs.append(
						(
							(interval_pairs[-1][0][1],interval_pairs[-1][0][1]+size),
							(interval_pairs[-1][1][1],interval_pairs[-1][1][1]+size),
						)
					)
				interval_pairs.append(
						(
							(interval_pairs[-1][0][1],interval_pairs[-1][0][1]+tSize-tEnd),
							(interval_pairs[-1][1][1],interval_pairs[-1][1][1]+qSize-qEnd),
						)
					)
				assert interval_pairs[-1][0][1]==tSize
				assert interval_pairs[-1][1][1]==qSize

			self.tmp_chain_sequences.append(
					ChainSequence(
							interval_pairs=interval_pairs,
							sampling_step=self._sampling_step,
							invert=self._invert,
							name1=tName,
							name2=qName,
						)
				)

		self._chain_sequences=collections.OrderedDict(
				[
					(chain_sequence.name1,chain_sequence) for chain_sequence in tmp_chain_sequences
				]
			)


	def zero_based_transl(self, chromosome, coordinate):
		self._chain_sequences[chromosome].zero_based_transl(coordinate)

	def one_based_transl(self, chromosome, coordinate):
		self._chain_sequences[chromosome].one_based_transl(coordinate)

	def get_fasta_index(self):
		#TODO: method for creating FAI object from a Chain object
		pass

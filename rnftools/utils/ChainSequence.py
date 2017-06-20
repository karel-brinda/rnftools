import sys
import traceback


class ChainSequence:
	def __init__(
			self,
			interval_pairs,
			sampling_step,
			invert=False,
			name1=None,
			name2=None,
	):

		self._sampling_step = sampling_step
		self._invert = invert
		self._name1 = name1
		self._name2 = name2

		self._interval_pairs = []

		last_left1, last_left2, last_right1, last_right2 = None, None, None, None

		for ((left1, right1), (left2, right2)) in interval_pairs:
			assert (last_right1 is None) or (last_right1 == left1)
			assert (last_right2 is None) or (last_right2 == left2)
			assert left1 <= right1
			assert left2 <= right2

			self._interval_pairs.append(((left1, right1), (left2, right2)))
			last_left1, last_right1, last_left2, last_right2 = left1, right1, left2, right2

		assert self._interval_pairs[0][0][0] == 0
		assert self._interval_pairs[0][1][0] == 0

		if self._invert:
			self._name1, self._name2 = self._name2, self._name1
			self._interval_pairs = [x[::-1] for x in self._interval_pairs]

		self._length1 = self._interval_pairs[-1][0][1]
		self._length2 = self._interval_pairs[-1][1][1]

		self._sampling = []
		p = 0
		for i in range(self._length1 // self._sampling_step + 1):
			coordinate = i * self._sampling_step
			while not ( \
							self._interval_pairs[p][0][0] <= coordinate and \
								coordinate < self._interval_pairs[p][0][1] \
					):
				p += 1
			self._sampling.append(p)

	@property
	def name1(self):
		return self._name1

	@property
	def name2(self):
		return self._name2

	@property
	def length1(self):
		return self._length1

	@property
	def length2(self):
		return self._length2

	def zero_based_transl(self, coordinate):
		assert isinstance(coordinate, int), coordinate
		assert 0 <= coordinate
		assert coordinate < self.length1
		i = coordinate // self._sampling_step
		p = self._sampling[i]
		while not (self._interval_pairs[p][0][0] <= coordinate and coordinate < self._interval_pairs[p][0][1]):
			p += 1
		(left1, right1), (left2, right2) = self._interval_pairs[p]
		assert right1 - left1 > 0
		slope = 1.0 * (right2 - left2) / (right1 - left1)
		translated_coordinate = int(round(left2 + slope * (coordinate - left1)))
		return translated_coordinate

	def one_based_transl(self, coordinate):
		assert 0 <= coordinate
		if coordinate == 0:
			return 0
		else:
			return self.zero_based_transl(coordinate - 1) + 1

import rnftools.mishmash
import rnftools.utils


class Sample:
	"""Class for a sample to be simulated.

	Args:
		name (str): Name of the sample.
		reads_in_tuple (int): Number of subreads in a read tuple.
		paired_end_mode (str): Mode for output files ("bwa" or "bfast").
	"""

	def __init__(
			self,
			name,
			reads_in_tuple,
			paired_end_mode="bwa",
	):
		self._name = name
		self._sources = []
		self._dir = name
		self._reads_in_tuple = reads_in_tuple
		self._pair_end_mode = paired_end_mode

		if self._reads_in_tuple == 1:
			self._fq_fns = [self._name + ".fq"]
			self._mode = "single-end"
		elif self._reads_in_tuple == 2 and self._pair_end_mode == "bwa":
			self._fq_fns = [self._name + ".1.fq", self._name + ".2.fq"]
			self._mode = "paired-end-bwa"
		elif self._reads_in_tuple == 2 and self._pair_end_mode == "bfast":
			self._fq_fns = [self._name + ".fq"]
			self._mode = "paired-end-bfast"

		rnftools.mishmash.add_sample(self)

		if paired_end_mode not in ["bwa", "bfast"]:
			rnftools.utils.error(
				"paired_end_mode must be 'bwa' or 'bfast'",
				program="RNFtools",
				subprogram="MIShmash",
				exception=ValueError,
			)

	def get_name(self):
		return self._name

	def get_dir(self):
		return self._dir

	def get_sources(self):
		return self._sources

	def get_input(self):
		return [source.fq_fn() for source in self._sources]

	def add_source(self, source):
		if self._reads_in_tuple != source.get_reads_in_tuple():
			message = "It is not possible to combine read tuples with different number of reads in a single sample. "
			"Details: name='{}', old ends='{}', new ends='{}', source='{}'.".format(
				self._name,
				self._reads_in_tuple,
				source.get_reads_in_tuple(),
				source,
			)

			rnftools.utils.error(
				message,
				program="RNFtools",
				subprogram="MIShmash",
				exception=ValueError,
			)

		self._sources.append(source)

	def clean(self):
		for x in self._fq_fns + [self._dir]:
			rnftools.utils.shell('rm -fR "{}"'.format(x))

	def fq_fns(self):
		return self._fq_fns

	def create_fq(self):
		fq_merger = rnftools.rnfformat.FqMerger(
			mode=self._mode,
			input_files_fn=[source.fq_fn() for source in self._sources],
			output_prefix=self._name,
		)
		fq_merger.run()

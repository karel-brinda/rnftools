import rnftools.rnfformat
import re

reg_lrn = re.compile(r"^([!-?A-^`-~]*)__([0-9a-f]+)__([!-?A-^`-~]+)__([!-?A-^`-~]*)$")
reg_prefix_part = re.compile(r"^[!-?A-^`-~]*$")
reg_id_part = re.compile(r"^[0-9a-f]+$")
reg_segmental_part = re.compile(r"^(?:(\([0-9FRN,]*\))(?:,(?!$)|$))+$")
reg_suffix_part = re.compile(r"^(?:((?:[a-zA-Z0-9]+:){0,1})\[([!-?A-Z\\^`-~]*)\](?:,(?!$)|$))+$")
reg_segment = re.compile(r"^\(([0-9]+),([0-9]+),([FRN]),([0-9]+),([0-9]+)\)$")
reg_comment = re.compile(r"^\[([!-?A-Z\\^`-~]*)\]$")
reg_extension = re.compile(r"^\[([!-?A-Z\\^`-~]*)\]$")


class Validator:
	"""Class for validation of RNF.

	Args:
		initial_read_tuple_name (str): Initial read tuple name to detect profile (widths).
		report_only_first (bool): Report only first occurrence of every error.
		warnings_as_errors (bool): Treat warnings as errors (error code).
	"""

	def __init__(
			self,
			initial_read_tuple_name,
			report_only_first=True,
			warnings_as_errors=False,
	):
		self.report_only_first = report_only_first
		self.reported_errors = set()
		self.error_has_been_reported = False
		self.warning_has_been_reported = False
		self.warnings_as_errors = warnings_as_errors

		self.rnf_profile = rnftools.rnfformat.RnfProfile(read_tuple_name=initial_read_tuple_name)

	def validate(self, read_tuple_name):
		"""Check RNF validity of a read tuple.

		Args:
			read_tuple_name (str): Read tuple name to be checked.s
		"""
		if reg_lrn.match(read_tuple_name) is None:
			self.report_error(
				read_tuple_name=read_tuple_name,
				error_name="wrong_read_tuple_name_structure",
				message="'{}' is not matched".format(reg_lrn),
			)
		else:
			parts = read_tuple_name.split("__")

			if reg_prefix_part.match(parts[0]) is None:
				self.report_error(
					read_tuple_name=read_tuple_name,
					error_name="wrong_prefix_part",
					message="'{}' is not matched".format(reg_prefix_part),
				)

			if reg_id_part.match(parts[1]) is None:
				self.report_error(
					read_tuple_name=read_tuple_name,
					error_name="wrong_id_part",
					message="'{}' is not matched".format(reg_id_part),
				)

			if reg_segmental_part.match(parts[2]) is None:
				self.report_error(
					read_tuple_name=read_tuple_name,
					error_name="wrong_segmental_part",
					message="'{}' is not matched".format(reg_segmental_part),
				)

			if reg_suffix_part.match(parts[3]) is None:
				self.report_error(
					read_tuple_name=read_tuple_name,
					error_name="wrong_suffix_part",
					message="'{}' is not matched".format(reg_suffix_part),
				)

			if not self.rnf_profile.check(read_tuple_name):
				self.report_error(
					read_tuple_name=read_tuple_name,
					error_name="wrong_profile",
					message="Read has a wrong profile (wrong widths). It should be: {} but it is: {}.".format(
						self.rnf_profile,
						rnftools.rnfformat.RnfProfile(read_tuple_name=read_tuple_name),
					),
					warning=True,
				)

	def get_return_code(self):
		"""Get final return code (0 = ok, 1=error appeared).
		"""
		if self.error_has_been_reported:
			return 1
		if self.warning_has_been_reported and self.warnings_as_errors:
			return 1

	def report_error(
			self,
			read_tuple_name,
			error_name,
			wrong="",
			message="",
			warning=False
	):
		"""Report an error.

		Args:
			read_tuple_name (): Name of the read tuple.
			error_name (): Name of the error.
			wrong (str): What is wrong. 
			message (str): Additional msessage to be printed.
			warning (bool): Warning (not an error).
		"""
		if (not self.report_only_first) or (error_name not in self.reported_errors):
			print("\t".join(["error" if warning == False else "warning", read_tuple_name, error_name, wrong, message]))
		self.reported_errors.add(error_name)
		if warning:
			self.warning_has_been_reported = True
		else:
			self.error_has_been_reported = True

import rnftools.rnfformat
import re

reg_lrn=re.compile(r"^([!-?A-^`-~]*)__([0-9a-f]+)__([!-?A-^`-~]+)__([!-?A-^`-~]*)$")
reg_prefix_part=re.compile(r"^[!-?A-^`-~]*$")
reg_id_part=re.compile(r"^[0-9a-f]+$")
reg_segmental_part=re.compile(r"^(?:(\([0-9FRN,]*\))(?:,(?!$)|$))+$")
reg_suffix_part=re.compile(r"^(?:((?:[a-zA-Z0-9]+:){0,1})\[([!-?A-Z\\^`-~]*)\](?:,(?!$)|$))+$")
reg_segment=re.compile(r"^\(([0-9]+),([0-9]+),([FRN]),([0-9]+),([0-9]+)\)$")
reg_comment=re.compile(r"^\[([!-?A-Z\\^`-~]*)\]$")
reg_extension=re.compile(r"^\[([!-?A-Z\\^`-~]*)\]$")


class Validator:
	def __init__(
				self,
				initial_read_name,
				report_only_first=True,
				warnings_as_errors=False,
			):
		self.report_only_first=report_only_first
		self.reported_errors=set()
		self.error_has_been_reported=False
		self.warning_has_been_reported=False
		self.warnings_as_errors=warnings_as_errors

		self.rnf_profile=rnftools.rnfformat.RnfProfile(read_name=initial_read_name)

	def validate(self,read_name):
		if reg_lrn.match(read_name) is None:
			self.report_error(
					read_name=read_name,
					error_name="wrong_read_name_structure",
					message="'{}' is not matched".format(reg_lrn),
				)
		else:
			parts=read_name.split("__")

			if reg_prefix_part.match(parts[0]) is None:
				self.report_error(
						read_name=read_name,
						error_name="wrong_prefix_part",
						message="'{}' is not matched".format(reg_prefix_part),
					)

			if reg_id_part.match(parts[1]) is None:
				self.report_error(
						read_name=read_name,
						error_name="wrong_id_part",
						message="'{}' is not matched".format(reg_id_part),
					)

			if reg_segmental_part.match(parts[2]) is None:
				self.report_error(
						read_name=read_name,
						error_name="wrong_segmental_part",
						message="'{}' is not matched".format(reg_segmental_part),
					)

			if reg_suffix_part.match(parts[3]) is None:
				self.report_error(
						read_name=read_name,
						error_name="wrong_suffix_part",
						message="'{}' is not matched".format(reg_suffix_part),
					)

			widths=self._widths_from_read_name(read_name)
			if not self.rnf_profile.check(read_name):
				self.report_error(
						read_name=read_name,
						error_name="wrong_profile",
						message="Read has a wrong profile (wrong widths). It should be: {} but it is: {}.".format(
								self.rnf_profile,
								rnftools.rnfformat.RnfProfile(read_name=read_name),
							),
						warning=True,
					)


	def get_return_code(self):
		if self.error_has_been_reported:
			return 1
		if self.warning_has_been_reported and self.warnings_as_errors:
			return 1

	def report_error(self,read_name,error_name,wrong="",message="",warning=False):
		if (not self.report_only_first) or (error_name not in self.reported_errors):
			print("\t".join(["error" if warning==False else "warning",read_name,error_name,wrong,message]))
		self.reported_errors.add(error_name)
		if warning:
			self.warning_has_been_reported=True
		else:
			self.error_has_been_reported=True

	def _widths_from_read_name(self,read_name):
		parts=read_name.split("__")
		segments=parts[2][1:-1].split("),(")
		int_widths=map(len,segments[0].split(","))
		return [len(parts[0]),len(parts[1])]+list(int_widths)

import re
import snakemake

from .Chain import *
from .FaIdx import *


def error(
		message,
		program=None,
		subprogram=None, exception=None,
):
	if exception != None:
		assert issubclass(exception, Exception)

	if program == None:
		program_part = ""
		subprogram_part = ""
	else:
		program_part = "[{}] ".format(program)
		if subprogram == None:
			subprogram_part = ""
		else:
			subprogram_part = "{}: ".format(subprogram)

	cprint(
		"".join([program_part, subprogram_part, "Error: ", message]),
		"red",
		attrs=['bold'],
	)

	if exception != None:
		raise exception(message)


def shell(
		cmd,
		remove_spaces=True,
		async=False,
		iterable=False,
		read=False,
):
	if remove_spaces:
		# print("removing spaces from command")
		cmd = re.sub(r'[ \t\f\v]+', ' ', cmd).strip()

	return snakemake.shell(
		cmd=cmd,
		async=async,
		iterable=iterable,
		read=read,
	)

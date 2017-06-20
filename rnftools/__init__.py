import builtins
import sys
import snakemake
import os
import re
from termcolor import colored, cprint

import rnftools.mishmash
import rnftools.lavender
import rnftools.rnfformat
import rnftools.utils

# version detection
try:
	import pkg_resources

	__version__ = pkg_resources.get_distribution("rnftools").version
except:
	__version__ = ""

DEFAULT_RNFTOOLS_CONF = {
	'print_info': sys.argv[0] == "snakemake",
}

# default configuration
try:
	RNFTOOLS_CONF = builtins.RNFTOOLS_CONF
except AttributeError:
	RNFTOOLS_CONF = {}

assert type(RNFTOOLS_CONF) is dict, "builtins.RNFTOOLS_CONF must be a dictionary"
for key in DEFAULT_RNFTOOLS_CONF.keys():
	if not key in RNFTOOLS_CONF:
		RNFTOOLS_CONF[key] = DEFAULT_RNFTOOLS_CONF[key]

# print info
if RNFTOOLS_CONF["print_info"]:
	message("", program="RNFtools")
	message("RNFtools", program="RNFtools")
	message("~~~~~~~~", program="RNFtools")
	message("Version:     {}".format(__version__), program="RNFtools")
	message("Web:         http://karel-brinda.github.io/rnftools/", program="RNFtools")
	message("Contact:     Karel Brinda, karel.brinda@univ-mlv.fr", program="RNFtools")
	message("Publication: K. Brinda, V. Boeva, G. Kucherov. RNF: a general framework to evaluate", program="RNFtools")
	message("             NGS read mappers, Bioinformatics 32(1), 2016 [DOI:10.1093/bioinformatics/btv524].",
		program="RNFtools")
	message("", program="RNFtools")


def include():
	return os.path.join(
		os.path.dirname(__file__),
		"include_all.snake",
	)


def input():
	return [
		rnftools.lavender.input(),
		rnftools.mishmash.input(),
	]


def message(message, program=None, subprogram=None):
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
		"".join([program_part, subprogram_part, message]),
		"blue",
		attrs=['bold'],
	)


def error(message, program=None, subprogram=None, exception=None):
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

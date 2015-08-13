import builtins
import sys

# switch off info for smbl
try:
	builtins.SMBL_CONF["print_info"]=False
except AttributeError:
	builtins.SMBL_CONF={"print_info":False}


import rnftools.mishmash
import rnftools.lavender
import rnftools.rnfformat
import rnftools.utils

import os

# version detection
try:
	import pkg_resources
	__version__=pkg_resources.get_distribution("rnftools").version
except:
	__version__=""

DEFAULT_RNFTOOLS_CONF = {
		'print_info': sys.argv[0]=="snakemake",
	}

# default configuration
try:
	RNFTOOLS_CONF = builtins.RNFTOOLS_CONF
except AttributeError:
	RNFTOOLS_CONF = {}

assert type(RNFTOOLS_CONF) is dict, "builtins.RNFTOOLS_CONF must be a dictionary"
for key in DEFAULT_RNFTOOLS_CONF.keys():
	if not key in RNFTOOLS_CONF:
		RNFTOOLS_CONF[key]=DEFAULT_RNFTOOLS_CONF[key]

# print info
if RNFTOOLS_CONF["print_info"]:
	import smbl.messages

	smbl.messages.message("",program="RNFtools")
	smbl.messages.message("RNFtools",program="RNFtools")
	smbl.messages.message("~~~~~~~~",program="RNFtools")
	smbl.messages.message("Version:     {}".format(__version__),program="RNFtools")
	smbl.messages.message("Web:         http://karel-brinda.github.io/rnftools/",program="RNFtools")
	smbl.messages.message("Contact:     Karel Brinda, karel.brinda@univ-mlv.fr",program="RNFtools")
	smbl.messages.message("Publication: K. Brinda, V. Boeva, G. Kucherov. RNF: a general framework to evaluate",program="RNFtools")
	smbl.messages.message("             NGS read mappers, arXiv:1504.00556 [q-bio.GN], 2015, accepted to Bioinformatics.",program="RNFtools")
	smbl.messages.message("",program="RNFtools")


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

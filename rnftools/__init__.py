import rnftools.mishmash
import rnftools.lavender
import rnftools.rnfformat
import smbl.messages
import os
import pkg_resources

version=pkg_resources.require("rnftools")[0].version

smbl.messages.message("",program="RNFtools")
smbl.messages.message("RNFtools",program="RNFtools")
smbl.messages.message("~~~~~~~~",program="RNFtools")
smbl.messages.message("Version:     {}".format(version),program="RNFtools")
smbl.messages.message("Web:         http://karel-brinda.github.io/rnftools/",program="RNFtools")
smbl.messages.message("Author:      Karel Brinda, karel.brinda@univ-mlv.fr",program="RNFtools")
smbl.messages.message("Publication: K.Brinda, V.Boeva, G.Kucherov, arXiv:1504.00556 ",program="RNFtools")
smbl.messages.message("                 [q-bio.GN], 2015, "+
	"accepted to Bioinformatics.",program="RNFtools")
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

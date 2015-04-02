import rnftools.mishmash
import rnftools.lavender
import rnftools.rnfformat
import smbl.messages
import os

smbl.messages.message("",program="RNFtools")
smbl.messages.message("Web:    http://github.com/karel-brinda/rnftools",program="RNFtools")
smbl.messages.message("Author: Karel Brinda, karel.brinda@gmail.com",program="RNFtools")
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

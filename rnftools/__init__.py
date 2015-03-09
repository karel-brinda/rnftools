import rnftools.mishmash
import rnftools.lavender
import rnftools.rnfformat

import os

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

__all__ = []


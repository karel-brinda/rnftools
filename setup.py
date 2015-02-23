import sys
from setuptools import setup , find_packages

if sys.version_info < (3,2):
	print("At least Python 3.2 is required.\n", file=sys.stderr)
	exit(1)

setup(
	name = 'MIShmash',
	packages = ['mishmash'],
	package_dir = {"mishmash" : "mishmash"},
	package_data = {"mishmash" : ["*.snake"] },
	version = '0.0.6',
	description = 'Program for simulating NGS read respecting the RNF read-name format',
	#long_description = """ \  """,
	install_requires=[
		'snakemake',
		'smbl',
		'pysam',
	],
	zip_safe=False,
	author = 'Karel Břinda',
	author_email = 'karel.brinda@gmail.com',
	url = 'https://github.com/karel-brinda/mishmash',
	scripts=['bin/rnf-join-fq.py'],
	license = "MIT",
	keywords = ['Snakemake', 'bioinformatics', 'read simulator', 'NGS'],
	classifiers = [
		"Development Status :: 3 - Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
)

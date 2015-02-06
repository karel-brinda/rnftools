import sys
from setuptools import setup , find_packages

if sys.version_info < (3,2):
	print("At least Python 3.2 is required.\n", file=sys.stderr)
	exit(1)

setup(
	name = 'LAVEnder',
	packages = ['lavender'],
	package_dir = {"lavender" : "lavender"},
	package_data = {"lavender" : ["*.snake"] },
	version = '0.0.1',
	description = 'Program for evaluating mappers of NGS reads keeping the RNF naming convention.',
	#long_description = """ \  """,
	install_requires=[
		'snakemake','smbl'
	],
	zip_safe=False,
	author = 'Karel Břinda',
	author_email = 'karel.brinda@gmail.com',
	url = 'https://github.com/karel-brinda/lavender',
	license = "MIT",
	#download_url = 'https://github.com/karel-brinda/mishmash/tarball/0.0.1',
	keywords = ['Snakemake', 'bioinformatics', 'mapper', 'evaluator'],
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

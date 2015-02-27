import sys

try:
	from setuptools import setup # , find_packages
except ImportError:
	print("Please install setuptools before installing rnftools.", file=sys.stderr)
	exit(1)

exec(open("rnftools/version.py").read())

setup(
	name = 'rnftools',
	packages = ['rnftools','rnftools.mishmash','rnftools.lavender'],
	package_dir = {
		"rnftools":"rnftools",
		"rnftools.mishmash":"rnftools/mishmash",
		"rnftools.lavender":"rnftools/lavender",
	},
	package_data = {
		"rnftools.mishmash" : ["*.snake"],
		"rnftools.lavender" : ["*.snake"],
	},
	version = __version__,
	description = 'RNF-compatible software for simulating NGS reads and evaluating mappers of NGS reads.',
	#long_description = """ \  """,
	install_requires=[
		'snakemake',
		'smbl',
		'pysam',
	],
        scripts=['bin/rnf-join-fq.py'],
	zip_safe=False,
	author = 'Karel Břinda',
	author_email = 'karel.brinda@gmail.com',
	url = 'http://github.com/karel-brinda/rnftools',
	license = "MIT",
	keywords = ['Snakemake', 'Bioinformatics', 'NGS Read Mapping', 'NGS Mappers Evaluation','Next-Generation Sequencing'],
	classifiers = [
		"Development Status :: 3 - Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 3",
		"Operating System :: Unix",
		"Topic :: Scientific/Engineering :: Bio-Informatics",
	],
)

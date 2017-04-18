#!/usr/bin/env python

import setuptools

import os
import sys

if sys.version_info < (3, 2):
	sys.exit('Minimum supported Python version is 3.2')

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()

# Get the current version
exec(open("rnftools/version.py").read())

setuptools.setup(
	name='RNFtools',

	version=VERSION,

	description='RNF framework for NGS: simulation of reads, evaluation of mappers, conversion of RNF-compliant data.',

	long_description=long_description,

	url='http://karel-brinda.github.io/rnftools/',

	author='Karel Brinda',
	author_email='karel.brinda@univ-mlv.fr',

	license='MIT',

	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: MIT License',
		'Natural Language :: English',
		'Programming Language :: Python :: 3 :: Only',
		'Operating System :: Unix',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
	],

	keywords='Bioinformatics, Computational Biology, Next-Generation Sequencing, Read Simulation, Mappers Evaluation',

	packages=[
		'rnftools',
		'rnftools.mishmash',
		'rnftools.lavender',
		'rnftools.rnfformat',
		'rnftools.utils',
	],

	install_requires=[
		'argparse',
		'beautifulsoup4',
		'pyfaidx',
		'pysam',
		'snakemake',
		'svg42pdf>=0.1.1',
		'termcolor',
	],

	package_data={
		'rnftools': [
			'*.snake',
		],
		'rnftools.mishmash': [
			'*.snake',
		],
		'rnftools.lavender': [
			'*.snake',
		],
	},

	entry_points={
		'console_scripts': [
			'rnftools = rnftools.scripts:rnftools_script',
		]
	},
)

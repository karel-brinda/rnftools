import re
import argparse
import os
import sys
import rnftools
import textwrap



# todo: examples of usages for every subcommand (using epilog)

################################################################
################################################################
##
## RNFTOOLS SUBCOMMANDS
##
################################################################
################################################################


def _add_shared_params(parser, unmapped_switcher=False):
	parser.add_argument(
		'-o', '--rnf-fastq',
		type=argparse.FileType('w+'),
		metavar='file',
		dest='fq_fo',
		required=True,
		help='Output FASTQ file (- for standard output).',
	)

	parser.add_argument(
		'-x', '--faidx',
		type=argparse.FileType('r'),
		metavar='file',
		dest='fai_fo',
		required=True,
		help="FAI index of the reference FASTA file (- for standard input). It can be created using 'samtools faidx'."
	)

	parser.add_argument(
		'-g', '--genome-id',
		type=int,
		metavar='int',
		dest='genome_id',
		default=1,
		help='Genome ID in RNF (default: 1).',
	)

	if unmapped_switcher:
		parser.add_argument(
			'-u', '--allow-unmapped',
			action='store_false',
			dest='allow_unmapped',
			help='Allow unmapped reads.',
		)


################################
# ART, MASON
################################
def sam2rnf(args):
	"""Convert SAM to RNF-based FASTQ with respect to argparse parameters.

	Args:
		args (...): Arguments parsed by argparse
	"""

	rnftools.mishmash.Source.recode_sam_reads(
		sam_fn=args.sam_fn,
		fastq_rnf_fo=args.fq_fo,
		fai_fo=args.fai_fo,
		genome_id=args.genome_id,
		number_of_read_tuples=10 ** 9,
		simulator_name=args.simulator_name,
		allow_unmapped=args.allow_unmapped,
	)


def add_sam2rnf_parser(subparsers, subcommand, help, description, simulator_name=None):
	"""Add another parser for a SAM2RNF-like command.

	Args:
		subparsers (subparsers): File name of the genome from which read tuples are created (FASTA file).
		simulator_name (str): Name of the simulator used in comments.
	"""

	parser_sam2rnf = subparsers.add_parser(subcommand, help=help, description=description)

	parser_sam2rnf.set_defaults(func=sam2rnf)

	parser_sam2rnf.add_argument(
		'-s', '--sam',
		type=str,
		metavar='file',
		dest='sam_fn',
		required=True,
		help='Input SAM/BAM with true (expected) alignments of the reads  (- for standard input).'
	)

	_add_shared_params(parser_sam2rnf, unmapped_switcher=True)

	parser_sam2rnf.add_argument(
		'-n', '--simulator-name',
		type=str,
		metavar='str',
		dest='simulator_name',
		default=simulator_name,
		help='Name of the simulator (for RNF).' if simulator_name is not None else argparse.SUPPRESS,
	)


################################
# CURESIM
################################

def curesim2rnf(args):
	rnftools.mishmash.CuReSim.recode_curesim_reads(
		rnf_fastq_fo=args.fq_fo,
		curesim_fastq_fo=args.curesim_fastq_fo,
		fai_fo=args.fai_fo,
		genome_id=args.genome_id,
		number_of_read_tuples=10 ** 9,
	)


def add_curesim_parser(subparsers, subcommand, help, description):
	parser_curesim2rnf = subparsers.add_parser(subcommand, help=help, description=description)

	parser_curesim2rnf.set_defaults(func=curesim2rnf)

	parser_curesim2rnf.add_argument(
		'-c', '--curesim-fastq',
		type=argparse.FileType('r'),
		metavar='file',
		dest='curesim_fastq_fo',
		required=True,
		help='CuReSim FASTQ file (- for standard input).',
	)

	_add_shared_params(parser_curesim2rnf, unmapped_switcher=False)


################################
# DWGSIM
################################

def dwgsim2rnf(args):
	rnftools.mishmash.DwgSim.recode_dwgsim_reads(
		dwgsim_prefix=args.dwgsim_prefix,
		fastq_rnf_fo=args.fq_fo,
		fai_fo=args.fai_fo,
		genome_id=args.genome_id,
		number_of_read_tuples=10 ** 9,
		# allow_unmapped=args.allow_unmapped,
		estimate_unknown_values=args.estimate_unknown_values,
	)


def add_dwgsim_parser(subparsers, subcommand, help, description):
	parser_dwgsim2rnf = subparsers.add_parser(subcommand, help=help, description=description)

	parser_dwgsim2rnf.set_defaults(func=dwgsim2rnf)

	parser_dwgsim2rnf.add_argument(
		'-p', '--dwgsim-prefix',
		type=str,
		metavar='str',
		required=True,
		dest='dwgsim_prefix',
		help='Prefix for DwgSim.',
	)

	parser_dwgsim2rnf.add_argument(
		'-e', '--estimate-unknown',
		action='store_true',
		dest='estimate_unknown_values',
		help='Estimate unknown values.',
	)

	# _add_shared_params(parser_dwgsim2rnf, unmapped_switcher=True)
	_add_shared_params(parser_dwgsim2rnf)


################################
# WGSIM
################################

def wgsim2rnf(args):
	rnftools.mishmash.WgSim.recode_wgsim_reads(
		rnf_fastq_fo=args.fq_fo,
		wgsim_fastq_1_fn=args.wgsim_fastq_1,
		wgsim_fastq_2_fn=args.wgsim_fastq_2,
		fai_fo=args.fai_fo,
		genome_id=args.genome_id,
		number_of_read_tuples=10 ** 9,
	)


def add_wgsim_parser(subparsers, subcommand, help, description):
	parser_wgsim2rnf = subparsers.add_parser(subcommand, help=help, description=description)

	parser_wgsim2rnf.set_defaults(func=wgsim2rnf)

	parser_wgsim2rnf.add_argument(
		'-1', '--wgsim-fastq-1',
		type=str,
		metavar='file',
		dest='wgsim_fastq_1',
		required=True,
		help='First WgSim FASTQ file (- for standard input)',
	)

	parser_wgsim2rnf.add_argument(
		'-2', '--wgsim-fastq-2',
		type=str,
		metavar='file',
		dest='wgsim_fastq_2',
		required=False,
		help='Second WgSim FASTQ file (in case of paired-end reads, default: none).',
		default=None,
	)
	_add_shared_params(parser_wgsim2rnf, unmapped_switcher=True)


################################
# MERGE
################################

def merge(args):
	mixer = rnftools.rnfformat.FqMerger(
		mode=args.m,
		input_files_fn=args.i,
		output_prefix=args.o,
	)

	mixer.run()


def add_merge_parser(subparsers, subcommand, help, description):
	parser_merge = subparsers.add_parser(
		subcommand,
		help=help,
		description=description,
		epilog=textwrap.dedent("""\
				Source RNF-FASTQ files should satisfy the following conditions:
					1) Each file contains only reads corresponding to one genome (with the
					   same genome id).
					2) All files contain reads of the same type (single-end / paired-end).
					3) Reads with more reads per tuple (e.g., paired-end) have '/1', etc.
					   in suffix (for identification of nb of read).
				"""
		),
		formatter_class=argparse.RawTextHelpFormatter,
	)
	parser_merge.set_defaults(func=merge)

	parser_merge.add_argument(
		'-i',
		required=True,
		metavar='inp',
		nargs='+',
		help='input FASTQ files',
	)

	parser_merge.add_argument(
		'-m',
		required=True,
		metavar='mode',
		choices=rnftools.rnfformat.ALLOWED_MODES,
		# type=lambda x: is_valid_mode(parser,x),
		help='mode for mergeing files (single-end / paired-end-bwa / paired-end-bfast)',
	)

	parser_merge.add_argument(
		'-o',
		metavar='out',
		required=True,
		help='output prefix',
	)


################################
# SAM=>ROC
################################

def sam2roc(args):
	assert args.allowed_delta >= 0
	cmd = """
				rnftools sam2es -d {allowed_delta} -i "{sam_fn}" -o -| \
				rnftools es2et -i - -o - | \
				rnftools et2roc -i - -o "{roc_fn}"
			""".format(
		allowed_delta=args.allowed_delta,
		sam_fn=args.sam_fn,
		roc_fn=args.roc_fn,
	)
	print("Called command: ", os.linesep, re.sub(r'[ \t\f\v]+', ' ', cmd).strip(), file=sys.stderr)
	rnftools.shell(cmd)


def add_sam2roc_parser(subparsers, subcommand, help, description):
	parser_sam2roc = subparsers.add_parser(subcommand, help=help, description=description)

	parser_sam2roc.set_defaults(func=sam2roc)

	parser_sam2roc.add_argument(
		'-i', '--sam',
		type=str,
		metavar='file',
		dest='sam_fn',
		required=True,
		help='SAM/BAM with aligned RNF reads(- for standard input).',
	)

	parser_sam2roc.add_argument(
		'-o', '--roc',
		type=str,
		metavar='file',
		dest='roc_fn',
		required=True,
		help='Output ROC file (- for standard output).',
		default=None,
	)

	parser_sam2roc.add_argument(
		'-d', '--allowed-delta',
		type=int,
		metavar='int',
		dest='allowed_delta',
		required=False,
		help='Tolerance of difference of coordinates between true (i.e., expected) alignment and real alignment (very important parameter!) (default: {}).'.format(
			rnftools.lavender.DEFAULT_ALLOWED_DELTA),
		default=rnftools.lavender.DEFAULT_ALLOWED_DELTA,
	)


################################
# ES
################################

def sam2es(args):
	assert args.allowed_delta >= 0
	rnftools.lavender.Bam.bam2es(
		bam_fn=args.sam_fn,
		es_fo=args.es_fo,
		allowed_delta=args.allowed_delta,
	)


def add_sam2es_parser(subparsers, subcommand, help, description):
	parser_sam2es = subparsers.add_parser(subcommand, help=help, description=description)

	parser_sam2es.set_defaults(func=sam2es)

	parser_sam2es.add_argument(
		'-i', '--sam',
		type=str,
		metavar='file',
		dest='sam_fn',
		required=True,
		help='SAM/BAM with aligned RNF reads(- for standard input).',
	)

	parser_sam2es.add_argument(
		'-o', '--es',
		type=argparse.FileType('w+'),
		metavar='file',
		dest='es_fo',
		required=True,
		help='Output ES file (evaluated segments, - for standard output).',
		default=None,
	)

	parser_sam2es.add_argument(
		'-d', '--allowed-delta',
		type=int,
		metavar='int',
		dest='allowed_delta',
		required=False,
		help='Tolerance of difference of coordinates between true (i.e., expected) alignment and real alignment (very important parameter!) (default: {}).'.format(
			rnftools.lavender.DEFAULT_ALLOWED_DELTA),
		default=rnftools.lavender.DEFAULT_ALLOWED_DELTA,
	)


################################
# ET
################################

def es2et(args):
	rnftools.lavender.Bam.es2et(
		es_fo=args.es_fo,
		et_fo=args.et_fo,
	)


def add_es2et_parser(subparsers, subcommand, help, description):
	parser_es2et = subparsers.add_parser(subcommand, help=help, description=description)

	parser_es2et.set_defaults(func=es2et)

	parser_es2et.add_argument(
		'-i', '--es',
		type=argparse.FileType('r'),
		metavar='file',
		dest='es_fo',
		required=True,
		help='Input ES file (evaluated segments, - for standard input).',
	)

	parser_es2et.add_argument(
		'-o', '--et',
		type=argparse.FileType('w+'),
		metavar='file',
		dest='et_fo',
		required=True,
		help='Output ET file (evaluated read tuples, - for standard output).',
		default=None,
	)


################################
# ROC
################################

def et2roc(args):
	rnftools.lavender.Bam.et2roc(
		et_fo=args.et_fo,
		roc_fo=args.roc_fo,
	)


def add_et2roc_parser(subparsers, subcommand, help, description):
	parser_et2roc = subparsers.add_parser(subcommand, help=help, description=description)

	parser_et2roc.set_defaults(func=et2roc)

	parser_et2roc.add_argument(
		'-i', '--et',
		type=argparse.FileType('r'),
		metavar='file',
		dest='et_fo',
		required=True,
		help='Input ET file (evaluated read tuples, - for standard input).',
	)

	parser_et2roc.add_argument(
		'-o', '--roc',
		type=argparse.FileType('w+'),
		metavar='file',
		dest='roc_fo',
		required=True,
		help='Output ROC file (evaluated reads, - for standard output).',
		default=None,
	)


################################
# CHECK
################################

def _max_version(versions):
	from distutils.version import LooseVersion as V

	max_ver = None
	for version in versions:
		version = str(version)
		if max_ver is None:
			max_ver = version
		else:
			if V(version) > V(max_ver):
				max_ver = version
	return max_ver


def _print_package_line(package, inst_ver, comp, av_ver):
	print('{:20} {:15} {:4} {:20}'.format(package, inst_ver, comp, av_ver))


def check(args):
	import xmlrpc.client
	import pkg_resources
	import socket
	from distutils.version import LooseVersion as V

	print("Checking the latest version. Please note that this is a very experimental subcommand.")
	print()

	_print_package_line("package", "installed", "", "available")
	print("=" * 70)
	for package in ["RNFtools"]:
		# 1) get version on PyPI
		ver_web = "N/A"
		try:
			pypi_client = xmlrpc.client.ServerProxy('https://pypi.python.org/pypi')
			available = pypi_client.package_releases(package)
			ver_web_max = _max_version(available)
			if ver_web_max is not None:
				ver_web = ver_web_max
		except socket.error:
			ver_web = "N/A"

		# 2) get installed version
		try:
			ver_installed = pkg_resources.require(package)[0].version
		except:
			ver_installed = "N/A"

		# 3) print summary
		versions_available = "N/A" not in [ver_installed, ver_web]
		update = False
		if versions_available:
			if V(ver_installed) > V(ver_web):
				comparison = ">"
			elif V(ver_installed) == V(ver_web):
				comparison = "="
			else:
				comparison = "<"
				update = True
		else:
			comparison = " "

		# pkg_info = '{dist.project_name} {dist.version}'.format(dist=dist)
		_print_package_line(package, ver_installed, comparison, ver_web)

	if update:
		print()
		print("An update is available. You can install the latest version using by")
		print("   pip install --upgrade rnftools")
		print()
		print("Note that pip may be available under a different name (pip-3, pip-3.4, etc.).")
		print("Root account might be required for this operation.")
	elif versions_available:
		print()
		print("Your installation is up-to-date.")


def add_check_parser(subparsers, subcommand, help, description):
	parser_check = subparsers.add_parser(subcommand, help=help, description=description)
	parser_check.set_defaults(func=check)


################################
# PUBLICATION
################################

def publication(args):
	print()
	print("--------------------------------------------------------------------------------------------")
	print("  K. Brinda, V. Boeva, G. Kucherov: RNF: a general framework to evaluate NGS read mappers.")
	print("        Bioinformatics (2016) 32(1): 136-139. [DOI:10.1093/bioinformatics/btv524].")
	print("--------------------------------------------------------------------------------------------")
	print()
	print("@article{rnftools,")
	print("\tauthor  = {B{\\v r}inda, Karel AND Boeva, Valentina AND Kucherov, Gregory},")
	print("\ttitle   = {RNF: a general framework to evaluate NGS read mappers},")
	print("\tjournal = {Bioinformatics},")
	print("\tyear    = {2016},")
	print("\tnumber  = {1},")
	print("\tvolume  = {32},")
	print("\tpmid    = {26353839},")
	print("\tdoi     = {10.1093/bioinformatics/btv524}")
	print("}")
	print()


def add_publication_parser(subparsers, subcommand, help, description):
	parser_publication = subparsers.add_parser(subcommand, help=help, description=description)
	parser_publication.set_defaults(func=publication)


################################
# VALIDATE
################################

def validate(args):
	for i, x in enumerate(args.fastq_fn):
		if i % 4 == 0:
			read_tuple_name = x.partition("/")[0][1:]
			if i == 0:
				validator = rnftools.rnfformat.Validator(
					initial_read_tuple_name=read_tuple_name,
					report_only_first=args.report_only_first,
					warnings_as_errors=args.warnings_as_errors,
				)
				validator.validate(read_tuple_name=read_tuple_name)
	sys.exit(validator.get_return_code())


def add_validate_parser(subparsers, subcommand, help, description):
	parser_validate = subparsers.add_parser(subcommand, help=help, description=description)

	parser_validate.add_argument(
		'-i', '--fastq',
		type=argparse.FileType('r'),
		metavar='file',
		dest='fastq_fn',
		required=True,
		help='FASTQ file to be validated.',
	)

	parser_validate.add_argument(
		'-w', '--warnings-as-errors',
		action='store_true',
		dest='warnings_as_errors',
		help='Treat warnings as errors.',
	)

	parser_validate.add_argument(
		'-a', '--all-occurrences',
		action='store_false',
		dest='report_only_first',
		help='Report all occurrences of warnings and errors.',
	)

	parser_validate.set_defaults(func=validate)


################################
# LIFTOVER
################################


def liftover(args):
	input_format = None
	output_format = None

	if args.input_fn.lower()[-4:] == ".sam":
		input_format = "sam"
	if args.input_fn.lower()[-4:] == "-":
		input_format = "sam"
	elif args.input_fn.lower()[-4:] == ".bam":
		input_format = "bam"
	elif args.input_fn.lower()[-4:] == ".fq":
		input_format = "fq"
	elif args.input_fn.lower()[-4:] == ".fastq":
		input_format = "fastq"

	if args.output_fn.lower()[-4:] == ".sam":
		output_format = "sam"
	if args.output_fn.lower()[-4:] == "-":
		output_format = "sam"
	elif args.output_fn.lower()[-4:] == ".bam":
		output_format = "bam"
	elif args.output_fn.lower()[-4:] == ".fq":
		output_format = "fq"
	elif args.output_fn.lower()[-4:] == ".fastq":
		output_format = "fastq"

	if args.input_format is not None:
		assert args.input_format.lower() in ["sam", "bam", "fastq", "fq"]
		if args.input_format.lower() == "sam":
			input_format = "sam"
		elif args.input_format.lower() == "bam":
			input_format = "sam"
		elif args.input_format.lower() in ["fastq", "fq"]:
			input_format = "fq"

	if args.output_format is not None:
		assert args.output_format.lower() in ["sam", "bam", "fastq", "fq"]
		if args.output_format.lower() == "sam":
			output_format = "sam"
		elif args.output_format.lower() == "bam":
			output_format = "sam"
		elif args.output_format.lower() in ["fastq", "fq"]:
			output_format = "fq"

	if input_format == "fq":
		assert output_format == "fq"
	if input_format in ["sam", "bam"]:
		assert output_format in ["sam", "bam"]

	assert input_format is not None
	assert output_format is not None

	rnf_lifter = rnftools.rnfformat.RnfLifter(
		chain_fn=args.chain_fn,
		fai_fn=args.fai_fn,
		invert=args.invert,
	)

	if input_format == "fq" and output_format == "fq":
		with open(args.input_fn) as fastq_in_fo:
			with open(args.output_fn, "w+") as fastq_out_fo:
				rnf_lifter.lift_fastq(
					fastq_in_fo=fastq_in_fo,
					fastq_out_fo=fastq_out_fo,
					genome_id=args.genome_id,
				)
	else:
		rnf_lifter.lift_sam(
			sam_in_fn=args.input_fn,
			sam_out_fn=args.output_fn,
			genome_id=args.genome_id,
		)


def add_liftover_parser(subparsers, subcommand, help, description):
	parser_liftover = subparsers.add_parser(subcommand, help=help, description=description)

	parser_liftover.add_argument(
		'-c', '--chain',
		type=str,
		metavar='file',
		dest='chain_fn',
		help='Chain liftover file for coordinates transformation. [no transformation]',
	)

	parser_liftover.add_argument(
		'-g', '--genome-id',
		type=str,
		metavar='int',
		dest='genome_id',
		help='ID of genome to be transformed.',
		required=True,
	)

	parser_liftover.add_argument(
		'-x', '--faidx',
		type=str,
		metavar='file',
		dest='fai_fn',
		help='Fasta index of the reference sequence. [extract from chain file]',
	)

	parser_liftover.add_argument(
		'--invert',
		action='store_true',
		dest='invert',
		help='Invert chain file (transformation in the other direction).',
	)

	parser_liftover.add_argument(
		'--input-format',
		type=str,
		metavar='str',
		dest='input_format',
		help='Input format (SAM/BAM/FASTQ). [autodetect]',
	)

	parser_liftover.add_argument(
		'--output-format',
		type=str,
		metavar='str',
		dest='output_format',
		help='Output format (SAM/BAM/FASTQ).  [autodetect]',
	)

	parser_liftover.add_argument(
		'input_fn',
		type=str,
		metavar='input',
		help='Input file to be transformed (- for standard input).',
	)

	parser_liftover.add_argument(
		'output_fn',
		type=str,
		metavar='output',
		help='Output file to be transformed (- for standard output).',
	)

	parser_liftover.set_defaults(func=liftover)


################################################################
################################################################
##
## RNFTOOLS SCRIPT
##
################################################################
################################################################

def rnftools_script():
	# create the top-level parser

	# !!! all command without parameters must be listed here!
	if len(sys.argv) == 1 or (len(sys.argv) == 2 and (sys.argv[1] != "publication") and (sys.argv[1] != "check")):
		sys.argv.append("-h")

	parser = argparse.ArgumentParser(
		prog=os.path.basename(sys.argv[0]),
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=textwrap.dedent("""\
			==================================================
			RNFtools -  http://karel-brinda.github.io/rnftools
			--------------------------------------------------
			version: {}
			contact: Karel Brinda (karel.brinda@univ-mlv.fr)
			==================================================
			""".format(rnftools.__version__),
		)
	)

	subparsers = parser.add_subparsers(
	)

	subparsers.add_parser("", help="", description="")

	#
	# rnftools check
	#
	add_check_parser(
		subparsers=subparsers,
		subcommand="check",
		help="Check for the latest version.",
		description="Check if RNFtools are up-to-date.",
	)

	#
	# rnftools publication
	#
	add_publication_parser(
		subparsers=subparsers,
		subcommand="publication",
		help="Print information about the associated publication.",
		description="Print information about the associated publication.",
	)

	#
	# rnftools validate
	#
	add_validate_parser(
		subparsers=subparsers,
		subcommand="validate",
		help="Validate RNF names in a FASTQ file.",
		description="Validate RNF names in a FASTQ file.",
	)

	#
	# rnftools liftover
	#
	add_liftover_parser(
		subparsers=subparsers,
		subcommand="liftover",
		help="Liftover genomic coordinates in RNF names.",
		description="Liftover genomic coordinates in RNF names in a SAM/BAM files or in a FASTQ file.",
	)

	subparsers.add_parser("", help="", description="")
	subparsers.add_parser("", help="---------------------[MIShmash]---------------------", description="")

	#
	# rnftools sam2rnf
	#
	add_sam2rnf_parser(
		subparsers=subparsers,
		subcommand="sam2rnf",
		help="Convert a SAM/BAM file to RNF-FASTQ.",
		description="Convert a SAM/BAM file to RNF-FASTQ.",
		simulator_name=None
	)

	#
	# rnftools art2rnf
	#
	add_sam2rnf_parser(
		subparsers=subparsers,
		subcommand="art2rnf",
		help="Convert output of Art to RNF-FASTQ.",
		description="""Convert an Art SAM file to RNF-FASTQ. Note that Art produces non-standard SAM files
				and manual editation might be necessary. In particular, when a FASTA file contains comments,
				Art left them in the sequence name. Comments must be removed in their corresponding @SQ headers in the SAM file,
				otherwise all reads are considered to be unmapped by this script.
			""",
		simulator_name="art",
	)

	#
	# rnftools curesim2rnf
	#
	add_curesim_parser(
		subparsers=subparsers,
		subcommand="curesim2rnf",
		help="Convert output of CuReSim to RNF-FASTQ.",
		description="Convert a CuReSim FASTQ file to RNF-FASTQ.",
	)

	#
	# rnftools dwgsim2rnf
	#
	add_dwgsim_parser(
		subparsers=subparsers,
		subcommand="dwgsim2rnf",
		help="Convert output of DwgSim to RNF-FASTQ.",
		description="Convert a DwgSim FASTQ file (dwgsim_prefix.bfast.fastq) to RNF-FASTQ. ",
	)

	#
	# rnftools mason2rnf
	#
	add_sam2rnf_parser(
		subparsers=subparsers,
		subcommand="mason2rnf",
		help="Convert output of Mason to RNF-FASTQ.",
		description="Convert a Mason SAM file to RNF-FASTQ.",
		simulator_name="mason",
	)

	#
	# rnftools wgsim2rnf
	#
	add_wgsim_parser(
		subparsers=subparsers,
		subcommand="wgsim2rnf",
		help="Convert output of WgSim to RNF-FASTQ.",
		description="Convert WgSim FASTQ files to RNF-FASTQ.",
	)

	#
	# rnftools merge
	#
	add_merge_parser(
		subparsers=subparsers,
		subcommand="merge",
		help="Merge RNF-FASTQ files and convert the output to 'traditional' FASTQ.",
		description="todo",
	)

	subparsers.add_parser("", help="", description="")
	subparsers.add_parser("", help="---------------------[LAVEnder]---------------------", description="")

	#
	# rnftools sam2es
	#
	add_sam2es_parser(
		subparsers=subparsers,
		subcommand="sam2es",
		help="Convert SAM/BAM with reads in RNF to ES (evaluation of segments).",
		description="todo",
	)

	#
	# rnftools es2et
	#
	add_es2et_parser(
		subparsers=subparsers,
		subcommand="es2et",
		help="Convert ES to ET (evaluation of read tuples).",
		description="todo",
	)

	#
	# rnftools et2roc
	#
	add_et2roc_parser(
		subparsers=subparsers,
		subcommand="et2roc",
		help="Convert ET to ROC (final statistics).",
		description="todo",
	)

	#
	# rnftools sam2roc
	#
	add_sam2roc_parser(
		subparsers=subparsers,
		subcommand="sam2roc",
		help="Previous three steps in a single command.",
		description="todo",
	)

	args = parser.parse_args()
	args.func(args)

import sys
import argparse
import os
import rnftools

################################
################################
##
## Auxiliary functions
##
################################
################################


def _add_shared_params(parser, unmapped_switcher=False):
	parser.add_argument(
			'-o','--rnf-fastq',
			type=argparse.FileType('w+'),
			metavar='file',
			dest='fq_fo',
			required=True,
			help='Output FASTQ file (- for standard output).',
		)
	
	parser.add_argument(
			'-i','--fasta-index',
			type=argparse.FileType('r'),
			metavar='file',
			dest='fai_fo',
			required=True,
			help="FAI index of the reference FASTA file (- for standard input). It can be created using 'samtools faidx'."
		)
	
	parser.add_argument(
			'-g','--genome-id',
			type=int,
			metavar='int',
			dest='genome_id',
			default=1,
			help='Genome ID in RNF (default: 1).',
		)
	if unmapped_switcher:
		parser.add_argument(
				'-u','--allow-unmapped',
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
		number_of_read_tuples=10**9,
		simulator_name=args.simulator_name,
		allow_unmapped=args.allow_unmapped,
	)

def add_sam2rnf_parser(subparsers,subcommand,help,description,simulator_name=None):
	"""Add another parser for a SAM2RNF-like command.

	Args:
		subparsers (subparsers): File name of the genome from which read tuples are created (FASTA file).
		simulator_name (str): Name of the simulator used in comments.
	"""

	parser_sam2rnf = subparsers.add_parser(subcommand,help=help,description=description)
	parser_sam2rnf.set_defaults(func=sam2rnf)
	parser_sam2rnf.add_argument(
			'-s','--sam',
			type=str,
			metavar='file',
			dest='sam_fn',
			required=True,
			help='Input SAM/BAM with true (expected) alignments of the reads  (- for standard input).'
		)
	
	_add_shared_params(parser_sam2rnf,unmapped_switcher=True)

	parser_sam2rnf.add_argument(
			'-n','--simulator-name',
			type=str,
			metavar='str',
			dest='simulator_name',
			default=simulator_name,
			help='Name of the simulator (for RNF).' if simulator_name is not None else argparse.SUPPRESS,
		)

################################
# DWGSIM
################################

def dwgsim2rnf(args):
	rnftools.mishmash.DwgSim.recode_dwgsim_reads(
		dwgsim_prefix=args.dwgsim_prefix,
		fastq_rnf_fo=args.fq_fo,
		fai_fo=args.fai_fo,
		genome_id=args.genome_id,
		number_of_read_tuples=10**9,
		allow_unmapped=args.allow_unmapped,
	)

def add_dwgsim_parser(subparsers,subcommand,help,description):
	parser_dwgsim2rnf = subparsers.add_parser(subcommand,help=help,description=description)
	parser_dwgsim2rnf.set_defaults(func=dwgsim2rnf)
	parser_dwgsim2rnf.add_argument(
			'-p','--dwgsim-prefix',
			type=str,
			metavar='str',
			dest='dwgsim_prefix',
			help='Prefix for DwgSim.',
		)
	_add_shared_params(parser_dwgsim2rnf,unmapped_switcher=True)


################################
# WGSIM
################################

def wgsim2rnf(args):
	rnftools.mishmash.WgSim.recode_wgsim_reads(
		rnf_fastq_fo=args.fq_fo,
		wgsim_fastq_1=args.wgsim_fastq_1,
		wgsim_fastq_2=args.wgsim_fastq_2,
		fai_fo=args.fai_fo,
		genome_id=args.genome_id,
		number_of_read_tuples=10**9,
	)

def add_wgsim_parser(subparsers,subcommand,help,description):
	parser_wgsim2rnf = subparsers.add_parser(subcommand, help=help,description=description)
	parser_wgsim2rnf.set_defaults(func=wgsim2rnf)
	parser_wgsim2rnf.add_argument(
			'-1','--wgsim-fastq-1',
			type=str,
			metavar='file',
			dest='wgsim_fastq_1',
			required=True,
			help='First WgSim FASTQ file (- for standard input)',
		)
	parser_wgsim2rnf.add_argument(
			'-2','--wgsim-fastq-2',
			type=str,
			metavar='file',
			dest='wgsim_fastq_2',
			required=False,
			help='Second WgSim FASTQ file (in case of paired-end reads, default: none).',
			default=None,
		)
	_add_shared_params(parser_wgsim2rnf,unmapped_switcher=True)


################################
# CURESIMs
################################

def curesim2rnf(args):
	rnftools.mishmash.CuReSim.recode_curesim_reads(
		rnf_fastq_fo=args.fq_fo,
		curesim_fastq_fo=args.curesim_fastq_fo,
		fai_fo=args.fai_fo,
		genome_id=args.genome_id,
		number_of_read_tuples=10**9,
	)

def add_curesim_parser(subparsers,subcommand,help,description):
	parser_curesim2rnf = subparsers.add_parser(subcommand,help=help,description=description)
	parser_curesim2rnf.set_defaults(func=curesim2rnf)
	parser_curesim2rnf.add_argument(
			'-c','--curesim-fastq',
			type=argparse.FileType('r'),
			metavar='file',
			dest='curesim_fastq_fo',
			required=True,
			help='CuReSim FASTQ file (- for standard input).',
		)
	_add_shared_params(parser_curesim2rnf,unmapped_switcher=False)

################################
################################
##
## RNFTOOLS SCRIPT
##
################################
################################

def default_func(args):
	#print(args)
	pass

def rnftools_script():
	# create the top-level parser
	parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]))
	parser.add_argument('--version', action='version', version=rnftools.__version__)
	parser.set_defaults(func=default_func)
	subparsers = parser.add_subparsers(help='sub-command help')

	#
	# rnftools sam2rnf
	#
	add_sam2rnf_parser(
			subparsers=subparsers,
			subcommand="sam2rnf",
			help="Convert SAM file to RNF FASTQ.",
			description="Convert SAM file to RNF FASTQ.",
			simulator_name=None
		)

	#
	# rnftools art2rnf
	#
	add_sam2rnf_parser(
			subparsers=subparsers,
			subcommand="art2rnf",
			help="Convert output of Art to RNF FASTQ.",
			description="""Convert an Art SAM file to RNF FASTQ. Note that Art produces non-standard SAM files
				and manual editation might be necessary. In particular, when a FASTA file contains comments,
				Art left them in the sequence name. They must be removed in @SQ headers in the SAM file,
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
			help="Convert output of CuReSim to RNF FASTQ.",
			description="Convert a CuReSim FASTQ file to RNF FASTQ.",
		)

	#
	# rnftools dwgsim2rnf
	#
	add_dwgsim_parser(
			subparsers=subparsers,
			subcommand="dwgsim2rnf",
			help="Convert output of DwgSim to RNF FASTQ.",
			description="Convert a DwgSim FASTQ file (dwgsim_prefix.bfast.fastq) to RNF FASTQ. ",
		)

	#
	# rnftools mason2rnf
	#
	add_sam2rnf_parser(
			subparsers=subparsers,
			subcommand="mason2rnf",
			help="Convert output of Mason to RNF FASTQ.",
			description="Convert a Mason SAM file to RNF FASTQ.",
			simulator_name="mason",
		)

	#
	# rnftools wgsim2rnf
	#
	add_wgsim_parser(
			subparsers=subparsers,
			subcommand="wgsim2rnf",
			help="Convert output of WgSim to RNF FASTQ.",
			description="Convert WgSim FASTQ files to RNF FASTQ.",
		)

	args = parser.parse_args()
	args.func(args)

	####

	####
	parser_sam2mis = subparsers.add_parser('sam2mis', help='b help')

	####
	parser_mis2mir = subparsers.add_parser('mis2mir', help='b help')

	####
	parser_mir2roc = subparsers.add_parser('mir2roc', help='b help')

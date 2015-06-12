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
			help='Output FASTQ file.',
		)
	
	parser.add_argument(
			'-i','--fasta-index',
			type=argparse.FileType('r'),
			metavar='file',
			dest='fai_fo',
			required=True,
			help='FAI index of the reference FASTA file.'
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

def add_sam2rnf_parser(subparsers,subcommand,help,simulator_name=None):
	"""Add another parser for a SAM2RNF-like command.

	Args:
		subparsers (subparsers): File name of the genome from which read tuples are created (FASTA file).
		simulator_name (str): Name of the simulator used in comments.
	"""

	parser_sam2rnf = subparsers.add_parser(subcommand, help=help)
	parser_sam2rnf.set_defaults(func=sam2rnf)
	parser_sam2rnf.add_argument(
			'-s','--sam',
			type=str,
			metavar='file',
			dest='sam_fn',
			required=True,
			help='Input SAM/BAM with true alignments of the reads.'
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

def add_dwgsim_parser(subparsers,subcommand):
	parser_dwgsim2rnf = subparsers.add_parser(subcommand, help=help)
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

def add_wgsim_parser(subparsers,subcommand):
	parser_wgsim2rnf = subparsers.add_parser(subcommand, help=help)
	parser_wgsim2rnf.set_defaults(func=wgsim2rnf)
	parser_wgsim2rnf.add_argument(
			'-1','--wgsim-fastq-1',
			type=str,
			metavar='str',
			dest='wgsim_fastq_1',
			required=True,
			help='',
		)
	parser_wgsim2rnf.add_argument(
			'-2','--wgsim-fastq-2',
			type=str,
			metavar='str',
			dest='wgsim_fastq_2',
			required=False,
			help='',
			default=None,
		)
	_add_shared_params(parser_wgsim2rnf,unmapped_switcher=True)


################################
# CURESIMs
################################


################################
################################
##
## RNFTOOLS SCRIPT
##
################################
################################

def rnftools_script():
	# create the top-level parser
	parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]))
	parser.add_argument('--version', action='version', version=rnftools.__version__)
	subparsers = parser.add_subparsers()

	#
	# rnftools sam2rnf
	#
	add_sam2rnf_parser(
			subparsers=subparsers,
			subcommand="sam2rnf",
			help="Help message",
			simulator_name=None
		)
	#
	# rnftools mason2rnf
	#
	add_sam2rnf_parser(
			subparsers=subparsers,
			subcommand="mason2rnf",
			help="Help message",
			simulator_name="mason",
		)
	#
	# rnftools art2rnf
	#
	add_sam2rnf_parser(
			subparsers=subparsers,
			subcommand="art2rnf",
			help="Help message",
			simulator_name="art",
		)


	add_dwgsim_parser(
			subparsers=subparsers,
			subcommand="dwgsim2rnf",
		)

	add_wgsim_parser(
			subparsers=subparsers,
			subcommand="wgsim2rnf",
		)

	args = parser.parse_args()
	args.func(args)



	####

	####
	parser_curesim2rnf = subparsers.add_parser('curesim2rnf', help='b help')

	####

	####

	####
	####

	####
	parser_sam2mis = subparsers.add_parser('sam2mis', help='b help')

	####
	parser_mis2mir = subparsers.add_parser('mis2mir', help='b help')

	####
	parser_mir2roc = subparsers.add_parser('mir2roc', help='b help')






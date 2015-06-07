import sys
import argparse
import os
import rnftools

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
			'-i','--sam',
			type=str,
			metavar='file',
			dest='sam_fn',
			required=True,
			help='Input SAM/BAM with true alignments of the reads.'
		)
	parser_sam2rnf.add_argument(
			'-o','--fastq',
			type=argparse.FileType('w+'),
			metavar='file',
			dest='fq_fo',
			required=True,
			help='Output FASTQ file.',
		)
	parser_sam2rnf.add_argument(
			'-f','--fai',
			type=argparse.FileType('r'),
			metavar='file',
			dest='fai_fo',
			required=True,
			help='FAI index of the reference FASTA file.'
		)
	parser_sam2rnf.add_argument(
			'-g','--genome-id',
			type=int,
			metavar='int',
			dest='genome_id',
			default=1,
			help='Genome ID in RNF (default: 1).',
		)
	parser_sam2rnf.add_argument(
			'-u','--allow-unmapped',
			action='store_false',
			dest='allow_unmapped',
			help='Allow unmapped reads.',
		)
	parser_sam2rnf.add_argument(
			'-s','--simulator-name',
			type=str,
			metavar='str',
			dest='simulator_name',
			default=simulator_name,
			help='Name of the simulator (for RNF).' if simulator_name is not None else argparse.SUPPRESS,
		)

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

	args = parser.parse_args()
	args.func(args)



	####
	parser_art2rnf = subparsers.add_parser('art2rnf', help='b help')

	####
	parser_curesim2rnf = subparsers.add_parser('curesim2rnf', help='b help')

	####
	parser_dwgsim2rnf = subparsers.add_parser('dwgsim2rnf', help='b help')

	####
	parser_mason2rnf = subparsers.add_parser('mason2rnf', help='a help')

	####
	####

	####
	parser_sam2mis = subparsers.add_parser('sam2mis', help='b help')

	####
	parser_mis2mir = subparsers.add_parser('mis2mir', help='b help')

	####
	parser_mir2roc = subparsers.add_parser('mir2roc', help='b help')






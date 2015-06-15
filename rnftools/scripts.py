import sys
import argparse
import os
import rnftools
import textwrap

# todo: examples of usages for every subcommand (using epilog)

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
# CURESIM
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
# SAM=ROC
################################

def sam2roc(args):
	rnftools.lavender.Bam.bam2es(
			bam_fn=args.sam_fn,
			roc_fo=args.roc_fo,
			allowed_delta=args.allowed_delta,
		)

def add_sam2roc_parser(subparsers,subcommand,help,description):
	parser_sam2roc = subparsers.add_parser(subcommand,help=help,description=description)
	parser_sam2roc.set_defaults(func=sam2roc)
	parser_sam2roc.add_argument(
			'-i','--sam',
			type=str,
			metavar='file',
			dest='sam_fn',
			required=True,
			help='SAM/BAM with aligned RNF reads(- for standard input).',
		)
	parser_sam2roc.add_argument(
			'-o','--roc',
			type=argparse.FileType('w+'),
			metavar='file',
			dest='es_fo',
			required=True,
			help='Output ROC file (- for standard output).',
			default=None,
		)
	parser_sam2roc.add_argument(
			'-d','--allowed-delta',
			type=int,
			metavar='int',
			dest='allowed_delta',
			required=False,
			help='Tolerance of difference of coordinates between true (i.e., expected) alignment and real alignment (very important parameter!) (default: {}).'.format(rnftools.lavender.DEFAULT_ALLOWED_DELTA),
			default=rnftools.lavender.DEFAULT_ALLOWED_DELTA,
		)

################################
# ES
################################

def sam2es(args):
	rnftools.lavender.Bam.bam2es(
			bam_fn=args.sam_fn,
			es_fo=args.es_fo,
			allowed_delta=args.allowed_delta,
		)

def add_sam2es_parser(subparsers,subcommand,help,description):
	parser_sam2es = subparsers.add_parser(subcommand,help=help,description=description)
	parser_sam2es.set_defaults(func=sam2es)
	parser_sam2es.add_argument(
			'-i','--sam',
			type=str,
			metavar='file',
			dest='sam_fn',
			required=True,
			help='SAM/BAM with aligned RNF reads(- for standard input).',
		)
	parser_sam2es.add_argument(
			'-o','--es',
			type=argparse.FileType('w+'),
			metavar='file',
			dest='es_fo',
			required=True,
			help='Output ES file (evaluated segments, - for standard output).',
			default=None,
		)
	parser_sam2es.add_argument(
			'-d','--allowed-delta',
			type=int,
			metavar='int',
			dest='allowed_delta',
			required=False,
			help='Tolerance of difference of coordinates between true (i.e., expected) alignment and real alignment (very important parameter!) (default: {}).'.format(rnftools.lavender.DEFAULT_ALLOWED_DELTA),
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

def add_es2et_parser(subparsers,subcommand,help,description):
	parser_es2et = subparsers.add_parser(subcommand, help=help,description=description)
	parser_es2et.set_defaults(func=es2et)
	parser_es2et.add_argument(
			'-i','--es',
			type=argparse.FileType('r'),
			metavar='file',
			dest='es_fo',
			required=True,
			help='Input ES file (evaluated segments, - for standard input).',
		)
	parser_es2et.add_argument(
			'-o','--et',
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

def add_et2roc_parser(subparsers,subcommand,help,description):
	parser_et2roc = subparsers.add_parser(subcommand, help=help,description=description)
	parser_et2roc.set_defaults(func=et2roc)
	parser_et2roc.add_argument(
			'-i','--et',
			type=argparse.FileType('r'),
			metavar='file',
			dest='et_fo',
			required=True,
			help='Input ET file (evaluated read tuples, - for standard input).',
		)
	parser_et2roc.add_argument(
			'-o','--roc',
			type=argparse.FileType('w+'),
			metavar='file',
			dest='roc_fo',
			required=True,
			help='Output ROC file (evaluated reads, - for standard output).',
			default=None,
		)

################################
# PUBLICATION
################################

def publication(args):
	print()
	print("-------------------------------------------------------------------------------------------")
	print("  K. Brinda, V. Boeva, G. Kucherov: RNF: a general framework to evaluate NGS read mappers.")
	print("             arXiv:1504.00556 [q-bio.GN], accepted to Bioinformatics, 2015.")
	print("-------------------------------------------------------------------------------------------")
	print()
	print("@article{rnf,")
	print("\tauthor  = {B{\\v r}inda, Karel AND Boeva, Valentina AND Kucherov, Gregory},")
	print("\ttitle   = {RNF: a general framework to evaluate NGS read mappers},")
	print("\tyear    = {2015},")
	print("\tvolume  = {abs/1504.00556},")
	print("\tee      = {http://arxiv.org/abs/1504.00556},")
	print("}")
	print()


def add_publication_parser(subparsers,subcommand,help,description):
	parser_curesim2rnf = subparsers.add_parser(subcommand,help=help,description=description)
	parser_curesim2rnf.set_defaults(func=publication)


################################
################################
##
## RNFTOOLS SCRIPT
##
################################
################################

def default_func(args):
	pass

def rnftools_script():
	# create the top-level parser

	if len(sys.argv)==1 or (len(sys.argv)==2 and (sys.argv[1]!="publication")):
		sys.argv.append("-h")

	parser = argparse.ArgumentParser(
			prog=os.path.basename(sys.argv[0]),
			formatter_class=argparse.RawDescriptionHelpFormatter,
			description=textwrap.dedent("""
					================================================
					RNFtools -  http://rnftools.rtfd.org
					------------------------------------
					version: {}
					upgrade: pip3 install --upgrade rnftools
					contact: Karel Brinda (karel.brinda@univ-mlv.fr)
					================================================
					""".format(rnftools.__version__),
				)
		)
	#parser.add_argument('--version', action='version', version=rnftools.__version__)
	parser.set_defaults(func=default_func)
	subparsers = parser.add_subparsers(
			help='----------------------------------------------------',
		)

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
	# rnftools publication
	#
	add_publication_parser(
			subparsers=subparsers,
			subcommand="publication",
			help="Print information about the associated publication.",
			description="",
		)

	#
	# rnftools sam2es
	#
	add_sam2es_parser(
			subparsers=subparsers,
			subcommand="sam2es",
			help="todo",
			description="todo",
		)

	#
	# rnftools es2et
	#
	add_es2et_parser(
			subparsers=subparsers,
			subcommand="es2et",
			help="todo",
			description="todo",
		)

	#
	# rnftools et2roc
	#
	add_et2roc_parser(
			subparsers=subparsers,
			subcommand="et2roc",
			help="todo",
			description="todo",
		)

	#
	# rnftools sam2roc
	#
	add_sam2roc_parser(
			subparsers=subparsers,
			subcommand="sam2roc",
			help="todo",
			description="todo",
		)

	args = parser.parse_args()
	args.func(args)

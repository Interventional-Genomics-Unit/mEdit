# == Native Modules ==
from argparse import ArgumentParser as argp
from argparse import RawTextHelpFormatter
import textwrap
# == Installed Modules ==
from importlib.metadata import version


def parse_arguments():
	# -> === Launch argparse parser === <-
	parser = argp(
		prog='mEdit',
		description=f'version {version("meditability")}',
		# epilog="mEdit is pretty cool, huh? :)",
		usage='%(prog)s [options]',
		formatter_class=RawTextHelpFormatter
	)

	programs = parser.add_subparsers(
		title="== mEdit Programs ==",
		description=textwrap.dedent('''
		mEdit can be operated through a list of different programs.'''),
		dest="program",
	)
	# === Db Setup ===
	dbset_parser = programs.add_parser(
		'db_set',
		help=textwrap.dedent('''
			Setup the necessary background data to run mEdit'''),
		formatter_class=RawTextHelpFormatter
	)
	ref_db_parse = dbset_parser.add_argument_group("== Reference Database Pre-Processing ==")
	ref_db_parse.add_argument('-d',
							  dest='db_path',
							  default='.',
							  help=textwrap.dedent('''
	                          Provide the path where the "mEdit_database" 
	                          directory will be created ahead of the analysis.
	                          Requires ~6.5GB in-disk storage 
	                          [default: ./mEdit_database]'''))
	ref_db_parse.add_argument('-l',
							  dest='latest_reference',
							  action='store_true',
							  help=textwrap.dedent('''
							  Request the latest human genome reference as part
							  of mEdit database unpacking. This is especially 
							  recommended when running predictions on private 
							  genome assemblies. [default: False]'''))
	ref_db_parse.add_argument('-c',
							  dest='custom_reference',
							  help=textwrap.dedent('''
							  Provide the path to a custom human reference genome 
							  in FASTA format. ***Chromosome annotation must follow a
							    ">chrN" format (case sensitive)'''))
	ref_db_parse.add_argument('-p',
							  dest='threads',
							  default='1',
							  help=textwrap.dedent('''
	                          Provide the number of cores for parallel decompression
	                          of mEdit databases.
	                          '''))

	# === Editors List ===
	list_parser = programs.add_parser(
		'list',
		help=textwrap.dedent('''
				Prints the current set of editors available on mEdit'''),
		formatter_class=RawTextHelpFormatter
	)
	editors_list = list_parser.add_argument_group("== Available Editors and BEs ==")
	editors_list.add_argument('-d',
							  dest='db_path',
							  default='.',
							  help=textwrap.dedent('''
	                          Provide the path where the "mEdit_database"
	                          directory was created ahead of the analysis
	                          using the "db_set" program.
	                          [default: ./mEdit_database]''')
							  )
	editors_list.add_argument('-e',
							  dest='editors',
							  action='store_true',
							  help=textwrap.dedent('''
	                          Provides the current list of available editors on mEdit
	                           '''))
	editors_list.add_argument('-b',
							  dest='base_editors',
							  action='store_true',
							  help=textwrap.dedent('''
	                          Provides the current list of available base editors on mEdit 
	                          '''))

	# === Guide Prediction Program ===
	fguides_parser = programs.add_parser(
		'guide_prediction',
		help=textwrap.dedent('''
			The core mEdit program finds potential guides for
			variants specified on the input by searching a diverse set of
			editors.'''),
		formatter_class=RawTextHelpFormatter
	)
	in_out = fguides_parser.add_argument_group("== Input/Output Options ==")
	in_out.add_argument(
		'-i',
		dest='query_input',
		required=True,
		help=textwrap.dedent('''
			Path to plain text file containing the query (or set of queries) 
			of variant(s) for mEdit analysis. Must be a single nucleotide 
			variation. See --qtype for formatting options.
			''')
	)
	in_out.add_argument(
		'-o',
		dest='output',
		default='medit_analysis',
		help=textwrap.dedent('''
			Path to root directory where mEdit outputs will be stored 
			[default: mEdit_analysis_<jobtag>/]''')
	)
	in_out.add_argument('-d',
						dest='db_path',
						default='.',
						help=textwrap.dedent('''
	                    Provide the path where the "mEdit_database" 
	                    directory was created ahead of the analysis 
	                    using the "db_set" program. 
	                    [default: ./mEdit_database]''')
						)
	in_out.add_argument('-j',
						dest='jobtag',
						help=textwrap.dedent('''
	                    Provide the tag associated with the current mEdit job.
	                    mEdit will generate a random jobtag by default''')
						)
	run_params = fguides_parser.add_argument_group("== mEdit Core Parameters ==")
	run_params.add_argument(
		'-m',
		dest='mode',
		default='standard',
		choices=['standard', 'private'],
		help=textwrap.dedent('''
			The MODE option determines how mEdit will run your job. 
			[default = "standard"]
			[1-] "standard": will find and process guides based on a 
			reference human genome assembly along with a diverse set of 
			pangenomes from HPRC.
			[2-] "private": requires a private VCF file that will be 
			 processed for guide prediction.''')
	)
	run_params.add_argument(
		'-g',
		dest='private_genome',
		default=None,
		help=textwrap.dedent('''
			Provide a gunzip compressed VCF file to run mEdit’s 
			private mode''')
	)
	run_params.add_argument(
		'--qtype',
		dest='qtype_request',
		default='hgvs',
		choices=['hgvs', 'coord'],
		help=textwrap.dedent('''
			Set the query type provided to mEdit. [default = "hgvs"]
			[1-] "hgvs": must at least contain the Refseq identifier 
			followed by “:” and the commonly used HGVS nomenclature. 
			Example: NM_000518.5:c.114G>A
			[2-] "coord": must contain hg38 coordinates followed by 
			(ALT>REF). Alleles must be the plus strand.
			Example: chr11:5226778C>T\n''')
	)
	run_params.add_argument(
		'--editor',
		dest='editor_request',
		default='clinical',
		choices=['clinical', 'custom', 'user defined list'],
		help=textwrap.dedent('''
			Delimits the set of editors to be used by mEdit. 
			[default = "clinical"]
			Use the "medit list" prompt to access the arrays of editors currently 
			supported in each category.
			[1-] "clinical": a short list of clinically relevant editors 
			that are either in pre-clinical or clinical trials.
			[2-] "custom": select guide search parameters. This requires a
			 separate input of parameters : ‘pam’, ‘pamISfirst’,’guidelen’
			[3-] "user defined list": - Comma-separated list of editors''')
	)
	run_params.add_argument(
		'--be',
		dest='be_request',
		choices=['off', 'default', 'custom', 'user defined list'],
		default='default',
		help=textwrap.dedent('''
			Add this flag to make mEdit process base-editors. 
			[default = off]''')
	)
	run_params.add_argument(
		'--dry',
		dest='dry_run',
		action='store_true',
		help=textwrap.dedent('''
			Perform a dry run of mEdit.'''))
	cluster_opt = fguides_parser.add_argument_group("== SLURM Options ==")
	cluster_opt.add_argument(
		'-p',
		dest='parallel_processes',
		help=textwrap.dedent('''
				Most processes in mEdit can be submitted to SLURM.
				When submitting mEdit jobs to SLURM, the user can specify
				the number of parallel processes that will be sent to the 
				server [default = 1]''')
	)
	cluster_opt.add_argument(
		'--ncores',
		dest='ncores',
		default=2,
		help=textwrap.dedent('''
			Specify the number of cores through which each parallel process 
			will be computed. [default = 2]''')
	)

	cluster_opt.add_argument(
		'--maxtime',
		dest='maxtime',
		default='1:00:00',
		help=textwrap.dedent('''
			Specify the maximum amount of time allowed for each parallel job.
			Format example: 2 hours -> "2:00:00" [default = 1 hour]''')
	)

	# === Off Target Effect Program ===
	casoff_parser = programs.add_parser(
		'offtargets',
		help=textwrap.dedent('''
			Predict off-target effect for the guides found'''),
		formatter_class=RawTextHelpFormatter
	)
	offtarget_params = casoff_parser.add_argument_group("== Off-Target Parameters ==")
	offtarget_params.add_argument(
		'-r', help='Reference Genome')

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments

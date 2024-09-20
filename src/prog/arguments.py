# == Native Modules ==
from argparse import ArgumentParser as argp
from argparse import RawTextHelpFormatter
import textwrap
# == Installed Modules ==
#from importlib.metadata import version


def parse_arguments():
	# -> === Launch argparse parser === <-
	parser = argp(
		prog='mEdit',
		#description=f'version {version("meditability")}',
		# epilog="mEdit is pretty cool, huh? :)",
		usage='%(prog)s ',
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
	                          Requires ~3.2GB in-disk storage 
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
	ref_db_parse.add_argument('-t',
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
	# editors_list.add_argument(('--editors'),
	# 						  dest='editors',
	# 						  action='store_true',
	# 						  help=textwrap.dedent('''
	#                           Provides the current list of available editors on mEdit
	#                            '''))
	# editors_list.add_argument('-b',
	# 						  dest='base_editors',
	# 						  action='store_true',
	# 						  help=textwrap.dedent('''
	#                           Provides the current list of available base editors on mEdit
	#                           '''))

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
			of variant(s) for mEdit analysis. See --qtype for formatting options.
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
	# TODO please change private to vcf or custom vcf, we dont want to insinuate private clinical info
	run_params.add_argument(
		'-m',
		dest='mode',
		default='standard',
		choices=['fast', 'standard', 'private'],
		help=textwrap.dedent('''
			The MODE option determines how mEdit will run your job. 
			[default = "standard"]
			[1-] "fast": will find and process guides based only on one 
			reference human genome.
			[2-] "standard": will find and process guides based on a 
			reference human genome assembly along with a diverse set of 
			pangenomes from HPRC.
			[3-] "private": requires a private VCF file that will be 
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
	# TODO please allow for gene and rsid. I started this for you, see commented out section
	run_params.add_argument(
		'--qtype',
		dest='qtype_request',
		default='hgvs',
		choices=['hgvs', 'coord'], #,'gene','rsid'],
		help=textwrap.dedent('''
			Set the query type provided to mEdit. [default = "hgvs"]
			[1-] "hgvs": must at least contain the Refseq identifier 
			followed by “:” and the commonly used HGVS nomenclature. 
			Example: NM_000518.5:c.114G>A, NM_000518.5(HBB):c.114G>A(p.Trp38Ter), NG_000007.3:g.70838G>A, NC_000011.10:g.5226778C>T
			[2-] "coord": must contain hg38 1-based coordinates followed by 
			(ALT>REF). Alleles must be the plus strand.
			Example: chr11:5226778C>T
			[3-] "gene": Provide the gene name to find all the clinvar variants listed in this gene
			Example: HBB
			[4-] "rsid": Provide the dbSNP rsID
			Example: rs33974936

			\n''')
	)
	run_params.add_argument(
		'--editor',
		dest='editor_request',
		default='clinical',
		choices=['clinical','custom', 'user defined list'],
		help=textwrap.dedent('''
			Delimits the set of editors to be used by mEdit. 
			[default = "clinical"]
			Use the "medit list" prompt to access the arrays of editors currently 
			supported in each category.
			[1-] "clinical": a short list of clinically relevant editors 
			that are either in pre-clinical or clinical trials.
			[2-]"user defined list" : - Comma-separated list of editors
			[3-] "custom": select guide search parameters. This requires a
			 separate input of parameters : ‘pam’, ‘pamISfirst’,’guidelen’ & 'dsb_pos'
			''')
	)
	run_params.add_argument(
		'--be',
		dest='be_request',
		choices=['off', 'default', 'custom', 'user defined list'],
		default='default',
		help=textwrap.dedent('''
		Add this flag to allow mEdit process base-editors. [default = off]
		[1-] “off”: disable base editor guides searching.
		[2-] “default”: use generic ABE and CBE with ‘NGG’ PAM and 4-8 base editing window
		[3-] “custom”: : select base editor search parameters. 
		This requires a separate input of parameters :
		 ‘be_pam’, ‘be_pamISfirst’,’be_guidelen’,’be_win’,’target_base’,’result_base’
		[4-]"user defined list": - Comma-separated list chosen from the “medit list” of base editors	
			''')
	)
	run_params.add_argument(
		'--cutdist',
		dest='cutdist',
		default='7',
		help=textwrap.dedent('''
				Max allowable window a variant start position can be from
				the editor cut site. This option not available for base editors. 
				[default = 7]''')
	)
	run_params.add_argument(
		'-gl',
		dest='guidelen',
		default='20',
		help=textwrap.dedent('''
		endonuclease spacer length for a custom editor. [default =20]
		Can ONLY be used for a ‘custom’ editor''')
	)
	run_params.add_argument(
		'-bgl',
		dest='be_guidelen',
		default='20',
		help=textwrap.dedent('''
			base editor spacer length for a custom editor. [default =20]
			Can ONLY be used for a ‘custom’ base editor''')
	)
	## TODO can this be an option that doesn't need an argument followed? like only use it if its TRUE

	run_params.add_argument(
		'--pamISfirst',
		dest='pamISfirst',
		default=False,
		help=textwrap.dedent('''
			Whether the PAM site is 5’ of target site [default = False]. 
			Can ONLY be used for a ‘custom’ editor''')
	)
	run_params.add_argument(
		'--be_pamISfirst',
		dest='be_pamISfirst',
		default=False,
		help=textwrap.dedent('''
			Whether the PAM site is 5’ of target site [default = False]. 
			Can ONLY be used for a ‘custom’ base editor''')

	)
	run_params.add_argument(
		'-dsb',
		dest='dsb_pos',
		default=-3,
		help=textwrap.dedent('''
			Double strand cut site relative to pam. This can be a single integer with a blunt end 
			endonuclease or 2 integers separated by a single comma when using an endonuclease that 
			produces staggered end cuts. for example spCas9 would be “-3” and Cas12 is “18,22”.
			[default = -3]
			Can ONLY be used for ‘custom’ endonuclease''')

	)
	run_params.add_argument(
		'-be_win',
		dest='be_win',
		default='4,8',
		help=textwrap.dedent('''
			Two positive integers separated by a comma that represent the base editing window. 
			The numbering begins at the 5’ most end.
			ex. CBE window is “4,8"
			Can ONLY be used for ‘custom’ base editor''')
	)

	run_params.add_argument(
		'-tarb',
		dest='target_base',
		default='A',
		help=textwrap.dedent('''
			A single  base that the custom base editor will target.
			ex. ABE target base is “A”
			Can ONLY be used for ‘custom’ base editor''')


	)
	run_params.add_argument(
		'-resb',
		dest='result_base',
		default='G',
		help=textwrap.dedent('''
				A single  base that the custom base editor will change the target base to.
				ex. ABE result base is “G”
				Can ONLY be used for ‘custom’ base editor''')

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
		'offtarget',
		help=textwrap.dedent('''
			Predict off-target effect for the guides found'''),
		formatter_class=RawTextHelpFormatter
	)
	offtarget_params = casoff_parser.add_argument_group("== Off-Target Parameters ==")
	offtarget_params.add_argument(
		'--dry',
		dest='dry_run',
		action='store_true',
		help=textwrap.dedent('''
				Perform a dry run of mEdit.'''))

	off_in_out = casoff_parser.add_argument_group("== Input/Output Options ==")
	off_in_out.add_argument(
		'-o',
		dest='output',
		default='medit_analysis',
		help=textwrap.dedent('''
		Path to root directory where mEdit guide_prediction
		 outputs were stored. "medit offtarget" can't 
		 operate if this path is incorrect. [default: mEdit_analysis_<jobtag>/]
		 ''')
	)
	off_in_out.add_argument('-d',
							dest='db_path',
							default='.',
							help=textwrap.dedent('''
		                    Provide the path where the "mEdit_database" 
		                    directory was created ahead of the analysis 
		                    using the "db_set" program. 
		                    [default: ./mEdit_database]''')
							)
	off_in_out.add_argument('-j',
							dest='jobtag',
							required=True,
							help=textwrap.dedent('''
							Provide the tag associated with the desired 
							"medit guide_prediction" job ID.
							"mEdit offtarget" will use the path from the
							 OUTPUT option to access this JOBTAG.
							''')
							)
	off_in_out.add_argument('--select_editors',
							dest='select_editors',
							default='',
							help=textwrap.dedent('''
								Provide a comma-separated list to select which 
								editors should be analyzed for offtarget effect.
								[default: all] 
								''')
							)

	off_cluster_opt = casoff_parser.add_argument_group("== SLURM Options ==")
	off_cluster_opt.add_argument(
		'-p',
		dest='parallel_processes',
		help=textwrap.dedent('''
					Most processes in mEdit can be submitted to SLURM.
					When submitting mEdit jobs to SLURM, the user can specify
					the number of parallel processes that will be sent to the 
					server [default = 1]''')
	)
	off_cluster_opt.add_argument(
		'--ncores',
		dest='ncores',
		default=2,
		help=textwrap.dedent('''
				Specify the number of cores through which each parallel process 
				will be computed. [default = 2]''')
	)

	off_cluster_opt.add_argument(
		'--maxtime',
		dest='maxtime',
		default='1:00:00',
		help=textwrap.dedent('''
				Specify the maximum amount of time allowed for each parallel job.
				Format example: 2 hours -> "2:00:00" [default = 1 hour]''')
	)
	# TODO: Finish the other options at the user interface

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments

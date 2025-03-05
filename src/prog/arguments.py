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
		usage='%(prog)s ',
		formatter_class=RawTextHelpFormatter
	)

	programs = parser.add_subparsers(
		title="== mEdit Programs ==",
		description=textwrap.dedent('''
		mEdit can be operated through a list of different programs.\n'''),
		dest="program",
	)
	# === Db Setup ===
	dbset_parser = programs.add_parser(
		'db_set',
		help=textwrap.dedent('''
			Setup the necessary background data to run mEdit\n'''),
		formatter_class=RawTextHelpFormatter
	)
	ref_db_parse = dbset_parser.add_argument_group("== Reference Database Pre-Processing ==")
	ref_db_parse.add_argument('-d',
							  dest='db_path',
							  default='.',
							  help=textwrap.dedent('''
	                          Provide the path where the "mEdit_database" 
	                          directory will be created ahead of the analysis.
	                          Requires ~7.5GB in-disk storage 
	                          [default: ./mEdit_database]\n'''))
	# ref_db_parse.add_argument('-l',
	# 						  dest='latest_reference',
	# 						  action='store_true',
	# 						  help=textwrap.dedent('''
	# 						  Request the latest human genome reference as part
	# 						  of mEdit database unpacking. This is especially
	# 						  recommended when running predictions on private
	# 						  genome assemblies. [default: False]\n'''))
	ref_db_parse.add_argument('-c',
							  dest='custom_reference',
							  help=textwrap.dedent('''
							  Provide the path to a custom human reference
							  genome  in FASTA format.
							  ***Chromosome annotation must follow a
							  ">chrN" format (case sensitive)\n'''))
	ref_db_parse.add_argument('-t',
							  dest='threads',
							  default='1',
							  help=textwrap.dedent('''
	                          Provide the number of cores for parallel 
	                          decompression of mEdit databases.
	                          \n'''))
	# ref_db_parse.add_argument('-m',
	# 						  dest='db_mode',
	# 						  default='mini',
	# 						  help=textwrap.dedent('''
	# 	                          Database download mode: Either "mini", or "full".
	# 	                          By default, "mini" will unpack ~5GB of data, allowing users
	# 	                          to carry out most of mEdit's functionalities, except for off-target
	# 	                          analysis utilizing HPRC's pangenomes. The "full" mode[default: "mini"].
	# 	                          \n'''))

	# === Editors List ===
	list_parser = programs.add_parser(
		'list',
		help=textwrap.dedent('''
				Prints the current set of editors available on mEdit\n'''),
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
	                          [default: ./mEdit_database]\n''')
							  )
	# editors_list.add_argument(('--editors'),
	# 						  dest='editors',
	# 						  action='store_true',
	# 						  help=textwrap.dedent('''
	#                           Provides the current list of available editors on mEdit
	#                            \n'''))
	# editors_list.add_argument('-b',
	# 						  dest='base_editors',
	# 						  action='store_true',
	# 						  help=textwrap.dedent('''
	#                           Provides the current list of available base editors on mEdit
	#                           \n'''))

	# === Guide Prediction Program ===

	fguides_parser = programs.add_parser(
		'guide_prediction',
		help=textwrap.dedent('''
			The core mEdit program finds potential guides for
			variants specified on the input by searching a diverse set of
			editors.\n'''),
		formatter_class=RawTextHelpFormatter
	)
	in_out = fguides_parser.add_argument_group("== Input/Output Options ==")
	in_out.add_argument(
		'-i',
		dest='query_input',
		required=True,
		help=textwrap.dedent('''
			Path to plain text file containing the query (or set of queries) 
			of variant(s) for mEdit analysis. 
			See --qtype for formatting options.
			\n''')
	)
	in_out.add_argument(
		'-o',
		dest='output',
		default='medit_analysis',
		help=textwrap.dedent('''
			Path to root directory where mEdit outputs will be stored 
			[default: mEdit_analysis_<jobtag>/]\n''')
	)
	in_out.add_argument('-d',
						dest='db_path',
						default='.',
						help=textwrap.dedent('''
	                    Provide the path where the "mEdit_database" 
	                    directory was created ahead of the analysis 
	                    using the "db_set" program. 
	                    [default: ./mEdit_database]\n''')
						)
	in_out.add_argument('-j',
						dest='jobtag',
						help=textwrap.dedent('''
	                    Provide the tag associated with the current mEdit job.
	                    mEdit will generate a random jobtag by default\n''')
						)
	run_params = fguides_parser.add_argument_group("== mEdit Core Parameters ==")
	run_params.add_argument(
		'-m',
		dest='mode',
		default='standard',
		choices=['fast', 'standard', 'vcf'],
		help=textwrap.dedent('''
			The MODE option determines how mEdit will run your job. 
			[default = "standard"]
			[1-] "fast": will find and process guides based only on one 
			reference human genome.
			[2-] "standard": will find and process guides based on a 
			reference human genome assembly along with a diverse set of 
			pangenomes from HPRC.
			[3-] "vcf": requires a custom VCF file that will be 
			processed for guide prediction.\n''')
	)
	run_params.add_argument(
		'-v',
		dest='custom_vcf',
		default=None,
		help=textwrap.dedent('''
			Provide a gunzip compressed VCF file to run mEdit’s 
			vcf mode\n''')
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
			Example: NM_000518.5:c.114G>A, 
			NM_000518.5(HBB):c.114G>A(p.Trp38Ter), 
			NG_000007.3:g.70838G>A, 
			NC_000011.10:g.5226778C>T
			[2-] "coord": must contain hg38 1-based coordinates followed by 
			(ALT>REF). Alleles must be the plus strand.
			Example: chr11:5226778C>T
			\n\n''')
	)
	run_params.add_argument(
		'--editor',
		dest='editor_request',
		default='clinical',
		help=textwrap.dedent('''
			Delimits the set of editors to be used by mEdit. 
			[default = "clinical"]
			[1-] "clinical": this value calls for a short list of clinically
			relevant editors that are either in pre-clinical or 
			clinical trials.
			[2-]<user defined list>: The user provides a comma-separated 
			list of editors. Use the 'medit list' command to access the current 
			set of available editors.
			[3-] "custom": Apply custom guide search parameters. 
			This requires a separate input of parameters: 
			‘--pam’, ‘--pamisfirst’,’--guidelen’, and 'dsb_pos'
			\n''')
	)
	run_params.add_argument(
		'--be',
		dest='be_request',
		default='default',
		help=textwrap.dedent('''
		Add this flag to allow mEdit process base-editors. [default = off]
		[1-] “off”: disable base editor guides searching.
		[2-] “default”: use generic ABE and CBE with ‘NGG’ PAM 
		and 4-8 base editing window
		[3-] “custom”: : select base editor search parameters. 
		This requires a separate input of parameters :
		‘--pam’, ‘--pamisfirst’,’--guidelen’,’--edit_win’,
		’--target_base’, and ’--result_base’
		[4-]<user defined list>: The user provides a comma-separated 
		list of base editors. Use 'medit list' to access the current set of
		available editors.
		\n''')
	)
	run_params.add_argument(
		'--cutdist',
		dest='cutdist',
		default='7',
		help=textwrap.dedent('''
					Max allowable window a variant start position can be from
					the editor cut site. This option not available for 
					base editors. 
					[default = 7]\n''')
	)
	run_params.add_argument(
		'--dry',
		dest='dry_run',
		action='store_true',
		help=textwrap.dedent('''
				Perform a dry run of mEdit.\n'''))

	custom_options = fguides_parser.add_argument_group("== Custom Editor Options ==")
	custom_options.add_argument(
		'--pam',
		dest='pam',
		default='XXX',
		help=textwrap.dedent('''
	        Specifies the PAM sequence to be used for custom guide or 
	        base editor searches. Required if "--editor custom" or
	        "--be custom" is used.
	    \n''')
	)
	custom_options.add_argument(
		'--guidelen',
		dest='guide_length',
		default=-1,
		help=textwrap.dedent('''
	        Specifies the length of the guide sequence for both 
	        custom endonuclease and base editor searches.
	        Required if "--editor custom" or "--be custom" is used.
	    \n''')
	)
	custom_options.add_argument(
		'--pamisfirst',
		dest='pam_is_first',
		action='store_true',
		help=textwrap.dedent('''
	        Indicates if the PAM is positioned before the guide for both
	        custom endonuclease and base editor searches. 
	        Required if "--editor custom" or "--be custom" is used.
	    \n''')
	)
	custom_options.add_argument(
		'--dsb_pos',
		dest='dsb_position',
		default=-10000,
		help=textwrap.dedent('''
	        Double strand cut site relative to pam. This can be a single
	        integer with a blunt end endonuclease or 2 integers separated 
	        by a single comma when using an endonuclease that produces 
	        staggered end cuts. For example spCas9 would be “-3” and 
	        Cas12 is “18,22”.
	        Required if "--editor custom" is used.
	    \n''')
	)
	custom_options.add_argument(
		'--edit_win',
		dest='editing_window',
		default=(0, 0),
		help=textwrap.dedent('''
	        Specifies the size of the editing window for custom base editor
	        search. Two positive integers separated by a comma that represent
	        the base editing window. The numbering begins at the 5’ most end.
	        For example: CBE window is “4,8"
	        Required if "--be custom" is used.
	    \n''')
	)
	custom_options.add_argument(
		'--target_base',
		dest='target_base',
		choices=['A', 'C', 'G', 'T'],
		default='X',
		help=textwrap.dedent('''
	        Specifies the target base for base editor modification.
	        For example: ABE target base is “A”
	        Required if "--be custom" is used.
	    \n''')
	)
	custom_options.add_argument(
		'--result_base',
		dest='result_base',
		choices=['A', 'C', 'G', 'T'],
		default='X',
		help=textwrap.dedent('''
	        Specifies the base that the target base will be converted
	        to for base editor searches. For example: ABE result base is “G”
	        Required if "--be custom" is used.
	    \n''')
	)

	cluster_opt = fguides_parser.add_argument_group("== SLURM Options ==")
	cluster_opt.add_argument(
		'--cluster',
		dest='cluster_request',
		action='store_true',
		help=textwrap.dedent('''
					Request job submission through SLURM  [default = None]\n''')
	)
	cluster_opt.add_argument(
		'-p',
		dest='parallel_processes',
		default=1,
		help=textwrap.dedent('''
				Most processes in mEdit can be submitted to SLURM.
					When submitting mEdit jobs to SLURM, the user can specify
					the number of parallel processes that will be sent to the 
					server. Otherwise, if applied to a local machine, this will 
					still parallelize some processes. [default = 1]\n''')
	)
	cluster_opt.add_argument(
		'--ncores',
		dest='ncores',
		default=1,
		help=textwrap.dedent('''
			Specify the number of cores through which each parallel process 
			will be computed. [default = 2]\n''')
	)

	cluster_opt.add_argument(
		'--maxtime',
		dest='maxtime',
		default='1:00:00',
		help=textwrap.dedent('''
			Specify the maximum amount of time allowed for each parallel job.
			Format example: 2 hours -> "2:00:00" [default = 1 hour]\n''')
	)

	# === Off Target Effect Program ===
	casoff_parser = programs.add_parser(
		'offtarget',
		help=textwrap.dedent('''
			Predict off-target effect for the guides found\n'''),
		formatter_class=RawTextHelpFormatter
	)
	offtarget_params = casoff_parser.add_argument_group("== Off-Target Parameters ==")
	offtarget_params.add_argument(
		'--dry',
		dest='dry_run',
		action='store_true',
		help=textwrap.dedent('''
				Perform a dry run of mEdit.\n'''))

	off_in_out = casoff_parser.add_argument_group("== Input/Output Options ==")
	off_in_out.add_argument(
		'-o',
		dest='output',
		default='medit_analysis',
		help=textwrap.dedent('''
		Path to root directory where mEdit guide_prediction
		outputs were stored. "medit offtarget" can't 
		operate if this path is incorrect. [default: mEdit_analysis_<jobtag>/]
		 \n''')
	)
	off_in_out.add_argument('-d',
							dest='db_path',
							default='.',
							help=textwrap.dedent('''
		                    Provide the path where the "mEdit_database" 
		                    directory was created ahead of the analysis 
		                    using the "db_set" program. 
		                    [default: ./mEdit_database]\n''')
							)
	off_in_out.add_argument('-j',
							dest='jobtag',
							required=True,
							help=textwrap.dedent('''
							Provide the tag associated with the desired 
							"medit guide_prediction" job ID.
							"mEdit offtarget" will use the path from the
							 OUTPUT option to access this JOBTAG.
							\n''')
							)
	off_in_out.add_argument('--select_editors',
							dest='select_editors',
							default='',
							help=textwrap.dedent('''
								Provide a comma-separated list to select which 
								editors should be analyzed for offtarget effect.
								[default: all] 
								\n''')
							)
	off_in_out.add_argument('--dna_bulge',
							dest='dna_bulge',
							default=0,
							help=textwrap.dedent('''
									Sets the number of insertions in the off-target.
									[default: 0] 
									\n''')
							)
	off_in_out.add_argument('--rna_bulge',
							dest='rna_bulge',
							default=0,
							help=textwrap.dedent('''
										Sets the number of deletions in the off-target.
										[default: 0] 
										\n''')
							)
	off_in_out.add_argument('--max_mismatch',
							dest='max_mismatch',
							default=3,
							help=textwrap.dedent('''
											Sets the maximum allowable number of mismatches.
											[default: 3] 
											\n''')
							)


	off_cluster_opt = casoff_parser.add_argument_group("== SLURM Options ==")
	off_cluster_opt.add_argument(
		'--cluster',
		dest='cluster_request',
		action='store_true',
		help=textwrap.dedent('''
						Request job submission through SLURM  [default = None]\n''')
	)
	off_cluster_opt.add_argument(
		'-p',
		dest='parallel_processes',
		default=1,
		help=textwrap.dedent('''
					Most processes in mEdit can be submitted to SLURM.
					When submitting mEdit jobs to SLURM, the user can specify
					the number of parallel processes that will be sent to the 
					server. Otherwise, if applied to a local machine, this will 
					still parallelize some processes. [default = 1]\n''')
	)
	off_cluster_opt.add_argument(
		'--ncores',
		dest='ncores',
		default=1,
		help=textwrap.dedent('''
				Specify the number of cores through which each parallel 
				process will be computed. [default = 2]\n''')
	)

	off_cluster_opt.add_argument(
		'--maxtime',
		dest='maxtime',
		default='1:00:00',
		help=textwrap.dedent('''
				Specify the maximum amount of time allowed for each parallel job.
				Format example: 2 hours -> "2:00:00" [default = 1 hour]\n''')
	)
	# TODO: Finish the other options at the user interface

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments

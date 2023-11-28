# == Native Modules ==
from argparse import ArgumentParser as argp
from datetime import datetime
import secrets
import string
import pytz
# == Installed Modules ==
import yaml
# == Project Modules ==
from programs.guide_prediction import run as guide_prediction
from programs.db_set import run as db_set


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='mEdit',
		description='',
		epilog="mEdit is pretty cool, huh? :)",
		usage='%(prog)s [options]'
	)
	programs = parser.add_subparsers(
		title="mEdit Programs",
		description="mEdit operates under a few different programs. Pick one!",
		dest="program",
        )
	# casoff = parser.add_subparsers(
	# 	title="Predict Off-Targets",
	# 	description="Predict off-target effects for the guides found",
	# 	dest="program",
	# )
	# === Db Setup ===
	dbset_parser = programs.add_parser(
		'db_set',
		help='Setup the necessary background data to run mEdit')
	ref_genome = dbset_parser.add_argument_group("== Reference Genome Pre-Processing ==")
	ref_genome.add_argument('--ref',
	                        dest='db_path',
	                        default='medit_database',
	                        help='Provide the path where mEdit background data should be'
	                             ' stored ahead of the analysis [default: mEdit_database_<jobtag>/')
	# === Guide Prediction Program ===
	fguides_parser = programs.add_parser(
		'guide_prediction',
		help='The core mEdit program finds potential guides for '
		     'variants specified on the input by searching a diverse set of editors'
	)
	in_out = fguides_parser.add_argument_group("== Input/Output Options ==")
	in_out.add_argument(
		'-i',
		dest='query_input',
		required=True,
		help='Path to plain text file containing the query (or set of queries) of variant(s) '
		     'for mEdit analysis'
	)
	in_out.add_argument(
		'-o',
		dest='output',
		default='mEdit_analysis',
		help='Path to root directory where mEdit output will be stored [default: mEdit_analysis_<jobtag>/]'
	)
	run_params = fguides_parser.add_argument_group("== mEdit Core Parameters ==")
	run_params.add_argument(
		'-m',
		dest='mode',
		default='standard',
		choices=['standard', 'private'],
		help='The MODE option determines how mEdit will run your job. '
	                         'On "standard" it will find and process guides based on a reference human genome assembly '
	                         'along with a diverse set of pangenomes from HPRC. '
	                         'On "private" it will require a private VCF file and use it to process and find guides. '
	                         '[default = "standard"]'
	)
	run_params.add_argument(
		'-g',
		dest='private_genome',
		help='Provide a gunzip compressed VCF file to run mEdit’s private mode'
	)
	run_params.add_argument(
		'--qtype',
		dest='qtype_request',
		default='hgvs',
		choices=['hgvs', 'coord'],
		help='Set the query type provided to mEdit. '
		     'Available types in the current version: "hgvs" or "coord" [default = "hgvs"]'
	)
	run_params.add_argument(
		'--editor',
		dest='editor_request',
		default='clinical',
		choices=['clinical', 'custom', 'user defined list'],
		help='Pick which set of editors will be used in the mEdit run. '
		     'Options: \n'
		     '"clinical" - EXPLAIN WHAT IS IN HERE; '
		     '"custom" - EXPLAIN WHAT IS IN HERE  '
		     '"user defined list" - EXPLAIN WHAT IS IN HERE [default = "clinical"]'
	)
	run_params.add_argument(
		'--be',
		dest='bemode_request',
		choices=['off', 'default', 'custom', 'user defined list'],
		help='Add this flag to make mEdit process base-editors [default = off]'
	)

	# === Off Target Effect Program ===
	casoff_parser = programs.add_parser(
		'offtargets',
		help='Predict off-target effects for the guides found'
	)
	offtarget_params = casoff_parser.add_argument_group("== Off-Target Parameters ==")
	offtarget_params.add_argument('-r', help='Reference Genome')

	cluster_opt = parser.add_argument_group("== SLURM Options ==")
	cluster_opt.add_argument(
		'--ncores',
		dest='cores',
		default=1,
		help='When submitting mEdit jobs to an HPC workload manager, '
		     'the user can specify the number of cores through which the '
		     'different variants will be parallelized [default = 1]'
	)
	cluster_opt.add_argument(
		'--maxtime',
		dest='maxtime',
		default='1:00:00',
		help='Specify the maximum amount of time allowed for each parallel job on mEdit run [default = 1:00:00]'
	)
	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def date_tag():
	# Create a random string of 20 characters
	random_str = ''.join(secrets.choice(string.ascii_letters + string.digits) for _ in range(10))

	# Set the timezone to PST
	pst = pytz.timezone('America/Los_Angeles')
	# Get the current date and time
	current_datetime = datetime.now(pst)
	# Format the date as a string with day, hour, minute, and second
	formatted_date = f"{current_datetime.strftime('%y%m%d%H%M%S%f')}_{random_str}"

	return formatted_date


def write_yaml_to_file(py_obj, filename):
	with open(f'{filename}.yaml', 'w',) as f:
		yaml.dump(py_obj, f, sort_keys=False)
	print('Written to file successfully')


def main():
	# === Call argument parsing function ===
	args = parse_arguments()
	# mEdit Program
	program = args.program

	# Assign jobtag and run mode to config
	jobtag = date_tag()

	if program == "guide_prediction":
		guide_prediction(args, jobtag)

	# == Database Parameters
	if program == "db_set":
		db_set(args, jobtag)


if __name__ == "__main__":
	main()

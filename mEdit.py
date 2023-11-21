# == Native Modules ==
from argparse import ArgumentParser as argp
from datetime import datetime
import secrets
import string
import pytz
# == Installed Modules ==
import yaml
# == Project Modules ==


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
	# === Guide Prediction Program ===
	fguides_parser = programs.add_parser(
		'find_guides',
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
		help='Path to root directory where mEdit output will be stored [default: mEdit_analysis/]'
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
		default='all',
		choices=['all'],
		help='Pick which set of editors will be used in the mEdit run. '
	                         'Options: \n'
	                         '"all" - EXPLAIN WHAT IS IN HERE; '
	                         '"other_options" - EXPLAIN WHAT IS IN HERE [default = "all"]'
	)
	run_params.add_argument(
		'--be',
		dest='bemode_request',
		action='store_true',
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
	# == Load Run Parameters values ==
	output_path = args.output
	run_mode = args.mode
	private_genome = args.private_genome
	qtype = args.qtype_request
	editor = args.editor_request
	beflag = args.bemode_request
	# == Load SLURM-related values ==
	ncores = args.ncores
	maxtime = args.maxtime

	# === Input Checks ===
	# Check run mode
	if run_mode == 'private':
		if not private_genome:
			raise "Please provide a VCF input file to run mEdit's private mode"
	# Process BEmode input
	bemode = 'on' if beflag else 'off'

	# === Load template configuration file ===
	with open("config/guide_pred.yaml", 'r') as config_handle:
		config_template = yaml.safe_load(config_handle)
	with open("config/cluster.yaml", 'r') as cluster_handle:
		cluster_template = yaml.safe_load(cluster_handle)

	# === Assign Variables to Configuration File ===
	# Assign jobtag and run mode to config
	jobtag = date_tag()
	config_template['run_name'] = f"{run_mode}_{jobtag}"
	config_template['processing_mode'] = run_mode
	config_template['output_directory'] = output_path
	# Assign run parameters to config
	config_template['qtype'] = qtype
	config_template['editor'] = editor
	config_template['BEmode'] = bemode
	# Assign cluster options
	cluster_template['__default__']['cores'] = ncores
	cluster_template['__default__']['time'] = maxtime

	# Set export paths for dynamic YAML files
	dynamic_config_path = f"config/config_{jobtag}.yaml"
	dynamic_cluster_path = f"config/cluster_{jobtag}.yaml"


if __name__ == "__main__":
	main()

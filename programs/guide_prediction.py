# == Native Modules ==
from os.path import abspath
# == Installed Modules ==
import yaml
# == Project Modules ==


def run(args, jobtag):
	# == Load Run Parameters values ==
	output_path = abspath(args.output)
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
	with open("config/medit_guide_pred.yaml", 'r') as config_handle:
		config_template = yaml.safe_load(config_handle)
	with open("config/medit_cluster.yaml", 'r') as cluster_handle:
		cluster_template = yaml.safe_load(cluster_handle)

	# === Assign Variables to Configuration File ===
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

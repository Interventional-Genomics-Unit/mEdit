# == Native Modules ==
import subprocess
from os.path import abspath
import sys
# == Installed Modules ==
import yaml
# == Project Modules ==
from programs.medit_lib import compress_file
from programs.medit_lib import write_yaml_to_file
from programs.medit_lib import set_export


def guide_prediction(args, jobtag):
	# == Load Run Parameters values ==
	root_dir = abspath(args.output)
	mode = args.mode
	private_genome = args.private_genome
	qtype = args.qtype_request
	editor = args.editor_request
	beflag = args.bemode_request
	# == Load SLURM-related values ==
	ncores = args.ncores
	maxtime = args.maxtime

	# === Set export paths tied to the SMK pipeline ===
	vcf_dir_path = f"{root_dir}/{mode}/source_vcfs"
	config_dir_path = f"{root_dir}/config"
	vcf_filename = f"{jobtag}.vcf"
	# === Set export paths for dynamic YAML files ===
	dynamic_config_path = f"{config_dir_path}/config_{jobtag}.yaml"
	dynamic_cluster_path = f"{config_dir_path}/cluster_{jobtag}.yaml"
	# => Create sub-folder to host private VCF
	set_export(vcf_dir_path)
	set_export(config_dir_path)

	print(f'Maybe {vcf_dir_path} was created')
	# Process BEmode input
	bemode = 'on' if beflag else 'off'

	# === Load template configuration file ===
	with open("config/medit_guide_pred.yaml", 'r') as config_handle:
		config_template = yaml.safe_load(config_handle)
	with open("config/medit_cluster.yaml", 'r') as cluster_handle:
		cluster_template = yaml.safe_load(cluster_handle)

	# === Assign Variables to Configuration File ===
	config_template['run_name'] = f"{mode}_{jobtag}"
	config_template['processing_mode'] = mode
	config_template['output_directory'] = root_dir
	# Assign run parameters to config
	config_template['qtype'] = qtype
	config_template['editor'] = editor
	config_template['BEmode'] = bemode
	# Assign cluster options
	cluster_template['__default__']['cores'] = ncores
	cluster_template['__default__']['time'] = maxtime

	# === Input Checks ===
	# Check run mode
	if mode == 'private':
		if not private_genome:
			print("Please provide a VCF input file to run mEdit's private mode")
			sys.exit(1)
		# ---- VCF ID adjustment for private vcf run ----
		# => Import VCF file prefix information to config file
		config_template["vcf_id"] = jobtag
		# => Create a copy of the VCF in the internal mEdit directory
		cmd_copy_vcf = f"cp {private_genome} {vcf_dir_path}/{vcf_filename}"
		subprocess.run(cmd_copy_vcf, shell=True)
		# => Check VCF file compression and compress if necessary
		compress_file(f"{vcf_dir_path}/{vcf_filename}")

	# === Write YAML configs to mEdit Root Directory ===
	write_yaml_to_file(config_template, dynamic_config_path)
	write_yaml_to_file(cluster_template, dynamic_cluster_path)

	# === Invoke SMK Pipelines ===
	snakemake_command = f"snakemake --snakefile vcf_processing.smk -j {ncores} --configfile {dynamic_config_path} --use-conda -n"

	# Execute the Snakemake command using subprocess
	try:
		print("Calling VCF Processing pipeline with the following command:")
		print(f"{snakemake_command}")
		subprocess.run(snakemake_command, shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(f"Error: {e}")

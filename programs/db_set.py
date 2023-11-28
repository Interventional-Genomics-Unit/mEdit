# == Native Modules ==

# == Installed Modules ==
import yaml
# == Project Modules ==


def run(args, jobtag):
	# === Load template configuration file ===
	with open("../config/medit_database.yaml", 'r') as config_handle:
		config_db_template = yaml.safe_load(config_handle)

	# === Load Database Path ===
	db_path = args.db_path
	db_path_tag = f"{db_path}_{jobtag}"
	# === Assign Variables to Configuration File ===
	# Assign jobtag and run mode to config
	config_db_template['fasta_root_path'] = f"{db_path_tag}/pkl"


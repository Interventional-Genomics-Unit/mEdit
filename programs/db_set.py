# == Native Modules ==
from os.path import abspath
import subprocess
# == Installed Modules ==
import yaml
# == Project Modules ==


def run(args, jobtag):
	# === Load template configuration file ===
	with open("config/medit_database.yaml", 'r') as config_handle:
		config_db_template = yaml.safe_load(config_handle)

	# === Load Database Path ===
	db_path = abspath(args.db_path)
	db_path_tag = f"{db_path}/medit_db-{jobtag}"
	# === Assign Variables to Configuration File ===
	# == Assign jobtag and Fasta root path ==
	fasta_root_path = f"{db_path_tag}/{config_db_template['fasta_root_path']}"
	config_db_template['fasta_root_path'] = fasta_root_path
	# == Parse the Processed Tables folder and its contents ==
	processed_tables = f"{db_path_tag}/{config_db_template['processed_tables']}"
	config_db_template["processed_tables"] = processed_tables
	config_db_template["simple_tables"] = f"{processed_tables}/{config_db_template['simple_tables']}"
	config_db_template["hgvs_lookup"] = f"{processed_tables}/{config_db_template['hgvs_lookup']}"
	config_db_template["clinvar_update"] = f"{processed_tables}/{config_db_template['clinvar_update']}"
	config_db_template["refseq_table"] = f"{processed_tables}/{config_db_template['refseq_table']}"

	# == Parse the Raw Tables folder and its contents ==
	raw_tables = f"{db_path_tag}/{config_db_template['raw_tables']}"
	config_db_template["raw_tables"] = raw_tables
	config_db_template["clinvar_summary"] = f"{raw_tables}/{config_db_template['clinvar_summary']}"
	config_db_template["hpa"] = f"{raw_tables}/{config_db_template['hpa']}"
	config_db_template["gencode"] = f"{raw_tables}/{config_db_template['gencode']}"

	print(config_db_template)
	# === Download Data ===
	# == SeqRecord Pickles
	print("Downloading Genomic References")
	cmd_aws = f"aws s3 cp --recursive s3://meditdb/pkl.gz {fasta_root_path}"
	# subprocess.run(cmd_aws, shell=True)
	# == Processed Tables and Raw Tables
	print("Downloading Pre-Processed Background Data Sets")
	cmd_aws = f"aws s3 cp --recursive s3://meditdb/processed_tables.gz {processed_tables}.gz"
	# subprocess.run(cmd_aws, shell=True)
	cmd_aws = f"aws s3 cp --recursive s3://meditdb/raw_tables.gz {raw_tables}.gz"
	# subprocess.run(cmd_aws, shell=True)

	print("Decompressing Database Genomic References")
	cmd_gz = f"tar zxf {fasta_root_path}.tar.gz"
	# subprocess.run(cmd_gz, shell=True)
	cmd_gz = f"tar zxf {processed_tables}.tar.gz"
	# subprocess.run(cmd_gz, shell=True)
	cmd_gz = f"tar zxf {raw_tables}.tar.gz"
	# subprocess.run(cmd_gz, shell=True)

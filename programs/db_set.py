# == Native Modules ==
from os.path import abspath
import subprocess
# == Installed Modules ==
import yaml
# == Project Modules ==
from programs.medit_lib import (launch_shell_cmd,
                                set_export,
                                write_yaml_to_file)


def dbset(args):
	# === Load template configuration file ===
	with open("config/medit_database.yaml", 'r') as config_handle:
		config_db_template = yaml.safe_load(config_handle)

	# === Load Database Path ===
	db_path_full = f"{abspath(args.db_path)}/medit_database"
	config_db_dir_path = f"{db_path_full}/config_db"

	vcf_dir_path = f"{db_path_full}/standard/source_vcfs"
	config_db_path = f"{config_db_dir_path}/config_db.yaml"

	set_export(vcf_dir_path)
	set_export(config_db_dir_path)
	# === Assign Variables to Configuration File ===
	#   == Parent Database Path
	config_db_template['db_path'] = f"{db_path_full}"
	#   == Assign jobtag and Fasta root path ==
	fasta_root_path = f"{db_path_full}/{config_db_template['fasta_root_path']}"
	config_db_template['fasta_root_path'] = fasta_root_path
	#   == Parse the Processed Tables folder and its contents ==
	processed_tables = f"{db_path_full}/{config_db_template['processed_tables']}"
	config_db_template["processed_tables"] = f"{processed_tables}"
	config_db_template["simple_tables"] = f"{processed_tables}/{config_db_template['simple_tables']}"
	config_db_template["hgvs_lookup"] = f"{processed_tables}/{config_db_template['hgvs_lookup']}"
	config_db_template["clinvar_update"] = f"{processed_tables}/{config_db_template['clinvar_update']}"
	config_db_template["refseq_table"] = f"{processed_tables}/{config_db_template['refseq_table']}"

	#   == Parse the Raw Tables folder and its contents ==
	raw_tables = f"{db_path_full}/{config_db_template['raw_tables']}"
	config_db_template["raw_tables"] = f"{raw_tables}"
	config_db_template["clinvar_summary"] = f"{raw_tables}/{config_db_template['clinvar_summary']}"
	config_db_template["hpa"] = f"{raw_tables}/{config_db_template['hpa']}"
	config_db_template["gencode"] = f"{raw_tables}/{config_db_template['gencode']}"

	# === Write YAML configs to mEdit Root Directory ===
	write_yaml_to_file(config_db_template, config_db_path)

	# === Download Data ===
	#   == SeqRecord Pickles
	print("Downloading Database of Genomic References")
	launch_shell_cmd(f"aws s3 cp s3://medit.db/genome_pkl.tar.gz {db_path_full}")
	#   == HPRC VCF files Setup
	launch_shell_cmd(f"aws s3 cp --recursive s3://medit.db/hprc/ {vcf_dir_path}")
	#   == Processed Tables and Raw Tables
	print("Downloading Pre-Processed Background Data Sets")
	launch_shell_cmd(f"aws s3 cp s3://medit.db/processed_tables.tar.gz {db_path_full}")
	launch_shell_cmd(f"aws s3 cp s3://medit.db/raw_tables.tar.gz {db_path_full}")
	print("Decompressing Databases")
	launch_shell_cmd(f"gzip -d {db_path_full}/*.gz")
	launch_shell_cmd(f"tar -xf {db_path_full}/raw_tables.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/raw_tables.tar")
	launch_shell_cmd(f"tar -xf {db_path_full}/processed_tables.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/processed_tables.tar")
	launch_shell_cmd(f"tar -xf {db_path_full}/genome_pkl.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/genome_pkl.tar")

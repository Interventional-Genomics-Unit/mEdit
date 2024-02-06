# == Native Modules ==
from os.path import abspath
# == Installed Modules ==
import yaml
# == Project Modules ==
from prog.medit_lib import (compress_file,
							download_s3_objects,
							pickle_chromosomes,
							project_file_path,
							launch_shell_cmd,
							set_export,
							write_yaml_to_file)


def dbset(args):
	# === Load template configuration file ===
	config_path = project_file_path("smk.config", "medit_database.yaml")
	with open(config_path, 'r') as config_handle:
		config_db_template = yaml.safe_load(config_handle)

	# === Load Database Path ===
	db_path_full = f"{abspath(args.db_path)}/medit_database"
	config_db_dir_path = f"{db_path_full}/config_db"
	# == Load args
	threads = args.threads
	latest_genome_download = args.latest_reference
	custom_reference = args.custom_reference

	vcf_dir_path = f"{db_path_full}/standard/source_vcfs"
	config_db_path = f"{config_db_dir_path}/config_db.yaml"

	set_export(vcf_dir_path)
	set_export(config_db_dir_path)
	# === Assign Variables to Configuration File ===
	#   == Parent Database Path
	config_db_template['meditdb_path'] = f"{db_path_full}"
	#   == Assign jobtag and Fasta root path ==
	fasta_root_path = f"{db_path_full}/{config_db_template['fasta_root_path']}"
	config_db_template['fasta_root_path'] = fasta_root_path
	#   == Assign Editor pickles path ==
	config_db_template["editors"] = f"{db_path_full}/{config_db_template['editors']}"
	config_db_template["base_editors"] = f"{db_path_full}/{config_db_template['base_editors']}"
	config_db_template["models_path"] = f"{db_path_full}/{config_db_template['models_path']}"
	#   == Parse the Processed Tables folder and its contents ==
	processed_tables = f"{db_path_full}/{config_db_template['processed_tables']}"
	config_db_template["processed_tables"] = f"{processed_tables}"
	# config_db_template["simple_tables"] = f"{processed_tables}/{config_db_template['simple_tables']}"
	# config_db_template["hgvs_lookup"] = f"{processed_tables}/{config_db_template['hgvs_lookup']}"
	# config_db_template["clinvar_update"] = f"{processed_tables}/{config_db_template['clinvar_update']}"
	config_db_template["refseq_table"] = f"{processed_tables}/{config_db_template['refseq_table']}"

	#   == Parse the Raw Tables folder and its contents ==
	# raw_tables = f"{db_path_full}/{config_db_template['raw_tables']}"
	# config_db_template["raw_tables"] = f"{raw_tables}"
	# config_db_template["clinvar_summary"] = f"{raw_tables}/{config_db_template['clinvar_summary']}"
	# config_db_template["hpa"] = f"{raw_tables}/{config_db_template['hpa']}"
	# config_db_template["gencode"] = f"{raw_tables}/{config_db_template['gencode']}"

	# === Download Data ===
	#   == SeqRecord Pickles
	print("# ---*--- Database of Genomic References ---*---")
	download_s3_objects("medit.db", "genome_pkl", fasta_root_path)
	# == Download the latest human reference genome by request
	if latest_genome_download:
		download_s3_objects("medit.db", "latest_genome_ref", fasta_root_path)
		local_latest_ref_path = f"{fasta_root_path}/latest_hg38.fa"
		pickle_chromosomes(local_latest_ref_path, fasta_root_path)
		launch_shell_cmd(f"bgzip -c -@ {threads} {local_latest_ref_path} > {local_latest_ref_path}.gz",
						 message="Compressing human reference genome")
		launch_shell_cmd(f"rm {local_latest_ref_path}",
						 message="Cleaning up unused files")
		config_db_template["latest_reference"] = "True"
	if custom_reference:
		local_custom_ref_path = f"{fasta_root_path}/custom_reference.fa"
		launch_shell_cmd(f"cp {custom_reference} {local_custom_ref_path}",
						 message="Setting up custom human reference genome")
		pickle_chromosomes(local_custom_ref_path, fasta_root_path)
		launch_shell_cmd(f"bgzip -c -@ {threads} {local_custom_ref_path} > {local_custom_ref_path}.gz",
						 message="Compressing human reference genome")
		launch_shell_cmd(f"rm {local_custom_ref_path}",
						 message="Cleaning up unused files")
		config_db_template["custom_reference"] = "True"

	# === Write YAML configs to mEdit Root Directory ===
	write_yaml_to_file(config_db_template, config_db_path)

	#   == HPRC VCF files Setup
	download_s3_objects("medit.db", "hprc", vcf_dir_path)

	#   == Processed Tables and Raw Tables
	print("# ---*--- Pre-Processed Background Data Sets ---*---")
	download_s3_objects("medit.db", "processed_tables.tar.gz", db_path_full)
	# download_s3_objects("medit.db", "raw_tables.tar.gz", db_path_full)
	download_s3_objects("medit.db", "pkl.tar.gz", db_path_full)

	#   == Decompress tar.gz files in the database
	print("# ---*--- Unpacking Background Data ---*---")
	launch_shell_cmd(f"gzip -d {db_path_full}/processed_tables.tar.gz", verbose=False,
					 check_exist=f"{db_path_full}/processed_tables")
	launch_shell_cmd(f"gzip -d {db_path_full}/pkl.tar.gz", verbose=False,
					 check_exist=f"{db_path_full}/pkl")
	# launch_shell_cmd(f"tar -xf {db_path_full}/raw_tables.tar --directory={db_path_full}/ && "
	#                  f"rm {db_path_full}/raw_tables.tar", message="Unpacking background tables")
	launch_shell_cmd(f"tar -xf {db_path_full}/processed_tables.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/processed_tables.tar",
					 message="Unpacking background tables", check_exist=f"{db_path_full}/processed_tables")
	launch_shell_cmd(f"tar -xf {db_path_full}/pkl.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/pkl.tar", message="Unpacking models", check_exist=f"{db_path_full}/pkl")
	launch_shell_cmd(f"gzip -d {config_db_template['refseq_table']}.gz", verbose=False,
					 check_exist=f"{config_db_template['refseq_table']}")

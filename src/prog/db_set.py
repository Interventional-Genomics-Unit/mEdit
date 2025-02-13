# == Native Modules ==
from os.path import abspath
# == Installed Modules ==
import yaml
# == Project Modules ==
from prog.medit_lib import (compress_file,
							download_gdrive_folder,
							download_s3_objects,
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
	# latest_genome_download = args.latest_reference
	custom_reference = args.custom_reference

	# === Assign internal variables  ===
	vcf_dir_path = f"{db_path_full}/standard/source_vcfs"
	config_db_path = f"{config_db_dir_path}/config_db.yaml"

	# === Allocate Files
	# === Assign Variables to Configuration File ===
	#   == Parent Database Path
	config_db_template['meditdb_path'] = f"{db_path_full}"
	#   == Assign jobtag and Fasta root path ==
	fasta_root_path = f"{db_path_full}/{config_db_template['fasta_root_path']}"
	config_db_template['fasta_root_path'] = fasta_root_path
	# #   == Bed Files path
	# config_db_template["bed_path"] = f"{db_path_full}/{config_db_template['bed_path']}"
	#   == GuideScan Indices path
	config_db_template["gscan_indices_path"] = f"{db_path_full}/{config_db_template['gscan_indices_path']}"
	#   == Assign Editor pickles path ==
	config_db_template["editors"] = f"{db_path_full}/{config_db_template['editors']}"
	config_db_template["base_editors"] = f"{db_path_full}/{config_db_template['base_editors']}"
	config_db_template["models_path"] = f"{db_path_full}/{config_db_template['models_path']}"
	#   == Parse the Processed Tables folder and its contents ==
	processed_tables = f"{db_path_full}/{config_db_template['processed_tables']}"
	config_db_template["processed_tables"] = f"{processed_tables}"
	config_db_template["refseq_table"] = f"{processed_tables}/{config_db_template['refseq_table']}"

	set_export(vcf_dir_path)
	set_export(config_db_dir_path)
	set_export(db_path_full)
	set_export(fasta_root_path)

	# === Download Data ===
	#   == SeqRecord Pickles
	print("# ---*--- Processing Database of Genomic References ---*---")
	if not custom_reference:
		download_s3_objects("medit.db", "genome_pkl", fasta_root_path)

	if custom_reference:
		local_custom_ref_path = f"{fasta_root_path}/custom_reference.fa"
		launch_shell_cmd(f"cp {custom_reference} {local_custom_ref_path}",
						 message="--> Setting up custom human reference genome")
		launch_shell_cmd(f"gzip {local_custom_ref_path}",
						 message="--> Compressing custom human reference genome")
		config_db_template["sequence_id"] = "custom_reference"

	# === Write YAML configs to mEdit Root Directory ===
	write_yaml_to_file(config_db_template, config_db_path)

	#   == HPRC VCF files Setup
	download_s3_objects("medit.db", "hprc", vcf_dir_path)

	#   == Bed and GuideScan indices Setup
	# download_s3_objects("medit.db", "bed_files.tar.gz", db_path_full)
	download_s3_objects("medit.db", "gscan_indices.tar.gz", db_path_full)

	#   == Processed Tables
	print("# ---*--- Downloading Pre-Processed Background Data Sets ---*---")
	download_s3_objects("medit.db", "processed_tables.tar.gz", db_path_full)
	download_s3_objects("medit.db", "pkl.tar.gz", db_path_full)

	#  == Decompress tar.gz files in the database ==> Uses parallel pigz when available
	print("# ---*--- Decompressing Background Data ---*---")
	# launch_shell_cmd(f"decompress -d {config_db_template['bed_path']}.tar.gz", verbose=False,
	# 				 check_exist=f"{config_db_template['bed_path']}.tar", message="Decompressing Bed files archive...")
	launch_shell_cmd(f"decompress -d {config_db_template['gscan_indices_path']}.tar.gz", verbose=False,
					 check_exist=f"{config_db_template['gscan_indices_path']}.tar", message="Decompressing Gscan Indices archive...")
	launch_shell_cmd(f"decompress -d {config_db_template['processed_tables']}.tar.gz", verbose=False,
					 check_exist=f"{config_db_template['processed_tables']}.tar", message="Decompressing Processed Tables archive...")
	launch_shell_cmd(f"decompress -d {db_path_full}/pkl.tar.gz", verbose=False,
					 check_exist=f"{db_path_full}/pkl.tar", message="Decompressing Models archive...")

	# launch_shell_cmd(f"tar -xf {db_path_full}/bed_files.tar --directory={db_path_full}/ && "
	# 				 f"rm {db_path_full}/bed_files.tar",
	# 				 check_exist=f"{config_db_template['bed_path']}", message="Unpacking Bed files...")
	launch_shell_cmd(f"tar -xf {config_db_template['gscan_indices_path']}.tar --directory={db_path_full}/ && "
					 f"rm {db_path_full}/gscan_indices.tar",
					 check_exist=f"{config_db_template['gscan_indices_path']}", message="Unpacking Guide Scan Index...")
	launch_shell_cmd(f"tar -xf {db_path_full}/processed_tables.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/processed_tables.tar",
					 check_exist=f"{db_path_full}/processed_tables", message="Unpacking Processed Tables...")
	launch_shell_cmd(f"tar -xf {db_path_full}/pkl.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/pkl.tar",
					 check_exist=f"{db_path_full}/pkl", message="Unpacking Models...")
	launch_shell_cmd(f"decompress -d {config_db_template['refseq_table']}.gz", verbose=False,
					 check_exist=f"{config_db_template['refseq_table']}")

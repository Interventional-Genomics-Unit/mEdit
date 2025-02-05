# == Native Modules ==
from os.path import abspath
# == Installed Modules ==
import yaml
# == Project Modules ==
from prog.medit_lib import (compress_file,
							download_gdrive_folder,
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
	# latest_genome_download = args.latest_reference
	# custom_reference = args.custom_reference

	# === Assign internal variables  ===
	vcf_dir_path = f"{db_path_full}/standard/source_vcfs"
	config_db_path = f"{config_db_dir_path}/config_db.yaml"

	set_export(vcf_dir_path)
	set_export(config_db_dir_path)
	set_export(db_path_full)

	# === Pull values from config variables ===
	standard_ref_prefix = config_db_template["sequence_id"]
	# cloud_storage_prefix = config_db_template["cloud_storage"]

	# # === Download GDrive files ===
	# print("# ---*--- Deploying mEdit Database ---*---")
	# download_gdrive_folder(cloud_storage_prefix, db_path_full)
	#
	# print("# ---*--- Decompressing Files ---*---")
	# #   == Decompress tar.gz files in the database ==> Uses parallel pigz when available
	# launch_shell_cmd(f"decompress -d {db_path_full}/*.gz", verbose=False)
	#
	# #   == Open TAR containers
	# #      -> Processed Tables
	# launch_shell_cmd(f"tar -xf {db_path_full}/processed_tables.tar --directory={db_path_full}/ && "
	# 				 f"rm {db_path_full}/processed_tables.tar",
	# 				 message="Unpacking background tables", check_exist=f"{db_path_full}/processed_tables")
	#
	# #      -> Pickled Models
	# launch_shell_cmd(f"tar -xf {db_path_full}/pkl.tar --directory={db_path_full}/ && "
	# 				 f"rm {db_path_full}/pkl.tar", message="Unpacking models", check_exist=f"{db_path_full}/pkl")
	#
	# #      -> Pangenomes
	# launch_shell_cmd(f"tar -xf {db_path_full}/hprc.tar --directory={db_path_full}/ && "
	# 				 f"rm {db_path_full}/hprc.tar", message="Unpacking pangenomes", check_exist=f"{db_path_full}/hprc")
	#
	# #      -> Bed Files
	# launch_shell_cmd(f"tar -xf {db_path_full}/bed_files.tar --directory={db_path_full}/ && "
	# 				 f"rm {db_path_full}/bed_files.tar", message="Unpacking Bed Files", check_exist=f"{db_path_full}/bed_files")
	#
	# #      -> Reference Genome
	# launch_shell_cmd(f"tar -xf {db_path_full}/genome_pkl.tar --directory={db_path_full}/ && "
	# 				 f"rm {db_path_full}/genome_pkl.tar", message="Unpacking Reference Genome",
	# 				 check_exist=f"{db_path_full}/genome_pkl")



	# === Allocate Files
	# launch_shell_cmd(f"mv {db_path_full}/hprc/* {vcf_dir_path} ")
	#
	# === Assign Variables to Configuration File ===
	#   == Parent Database Path
	config_db_template['meditdb_path'] = f"{db_path_full}"
	#   == Assign jobtag and Fasta root path ==
	fasta_root_path = f"{db_path_full}/{config_db_template['fasta_root_path']}"
	config_db_template['fasta_root_path'] = fasta_root_path
	#   == Bed Files path
	config_db_template["bed_path"] = f"{db_path_full}/{config_db_template['bed_path']}"
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

	# === Download Data ===
	#   == SeqRecord Pickles
	print("# ---*--- Processing Database of Genomic References ---*---")
	skip_genome_pkl = download_s3_objects("medit.db", "genome_pkl", fasta_root_path)

	standard_ref_path = f"{fasta_root_path}/{standard_ref_prefix}.fa"
	if not skip_genome_pkl:
		launch_shell_cmd(f"bgzip -df -@ {threads} {standard_ref_path}.gz > {standard_ref_path}",
						 message="Decompressing human reference genome")
		pickle_chromosomes(standard_ref_path, fasta_root_path)
		launch_shell_cmd(f"bgzip -cf -@ {threads} {standard_ref_path} > {standard_ref_path}.gz")
		launch_shell_cmd(f"rm {standard_ref_path}",
						 message="Cleaning up unused files")
	# # == Download the latest human reference genome by request
	# if latest_genome_download:
	# 	download_s3_objects("medit.db", "latest_genome_ref", fasta_root_path)
	# 	local_latest_ref_path = f"{fasta_root_path}/latest_hg38.fa"
	# 	pickle_chromosomes(local_latest_ref_path, fasta_root_path)
	# 	launch_shell_cmd(f"bgzip -cf -@ {threads} {local_latest_ref_path} > {local_latest_ref_path}.gz",
	# 					 message="Compressing human reference genome")
	# 	launch_shell_cmd(f"rm {local_latest_ref_path}",
	# 					 message="Cleaning up unused files")
	# 	config_db_template["latest_reference"] = "True"
	# if custom_reference:
	# 	local_custom_ref_path = f"{fasta_root_path}/custom_reference.fa"
	# 	launch_shell_cmd(f"cp {custom_reference} {local_custom_ref_path}",
	# 					 message="Setting up custom human reference genome")
	# 	pickle_chromosomes(local_custom_ref_path, fasta_root_path)
	# 	launch_shell_cmd(f"bgzip -c -@ {threads} {local_custom_ref_path} > {local_custom_ref_path}.gz",
	# 					 message="Compressing human reference genome")
	# 	launch_shell_cmd(f"rm {local_custom_ref_path}",
	# 					 message="Cleaning up unused files")
	# 	config_db_template["custom_reference"] = "True"

	# === Write YAML configs to mEdit Root Directory ===
	write_yaml_to_file(config_db_template, config_db_path)

	#   == HPRC VCF files Setup
	download_s3_objects("medit.db", "hprc", vcf_dir_path)

	#   == Bed and GuideScan indices Setup
	download_s3_objects("medit.db", "bed_files.tar.gz", db_path_full)
	download_s3_objects("medit.db", "gscan_indices.tar.gz", db_path_full)

	#   == Processed Tables
	print("# ---*--- Downloading Pre-Processed Background Data Sets ---*---")
	download_s3_objects("medit.db", "processed_tables.tar.gz", db_path_full)
	download_s3_objects("medit.db", "pkl.tar.gz", db_path_full)

	#  == Decompress tar.gz files in the database ==> Uses parallel pigz when available
	print("# ---*--- Decompressing Background Data ---*---")
	launch_shell_cmd(f"decompress -d {vcf_dir_path}/*.gz", verbose=False,
					 check_exist=f"{vcf_dir_path}", message="Decompressing HPRC data archive...")
	launch_shell_cmd(f"decompress -d {config_db_template['bed_path']}.tar.gz", verbose=False,
					 check_exist=f"{config_db_template['bed_path']}.tar", message="Decompressing Bed files archive...")
	launch_shell_cmd(f"decompress -d {config_db_template['gscan_indices_path']}.tar.gz", verbose=False,
					 check_exist=f"{config_db_template['gscan_indices_path']}.tar", message="Decompressing Gscan Indices archive...")
	launch_shell_cmd(f"decompress -d {config_db_template['processed_tables']}.tar.gz", verbose=False,
					 check_exist=f"{config_db_template['processed_tables']}.tar", message="Decompressing Processed Tables archive...")
	launch_shell_cmd(f"decompress -d {db_path_full}/pkl.tar.gz", verbose=False,
					 check_exist=f"{db_path_full}/pkl.tar", message="Decompressing Models archive...")

	launch_shell_cmd(f"tar -xf {db_path_full}/gscan_indices.tar --directory={db_path_full}/ && "
					 f"rm {db_path_full}/gscan_indices.tar",
					 check_exist=f"{db_path_full}/gscan_indices", message="Unpacking Guide Scan Index...")
	launch_shell_cmd(f"tar -xf {db_path_full}/processed_tables.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/processed_tables.tar",
					 check_exist=f"{db_path_full}/processed_tables", message="Unpacking Processed Tables...")
	launch_shell_cmd(f"tar -xf {db_path_full}/pkl.tar --directory={db_path_full}/ && "
	                 f"rm {db_path_full}/pkl.tar",
					 check_exist=f"{db_path_full}/pkl", message="Unpacking Models...")
	launch_shell_cmd(f"tar -xf {db_path_full}/bed_files.tar --directory={db_path_full}/ && "
					 f"rm {db_path_full}/bed_files.tar",
					 check_exist=f"{db_path_full}/bed_files",message="Unpacking Bed files...")
	launch_shell_cmd(f"decompress -d {config_db_template['refseq_table']}.gz", verbose=False,
					 check_exist=f"{config_db_template['refseq_table']}")

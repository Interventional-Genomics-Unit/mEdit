# == Native Modules ==
import subprocess
import gzip
import shutil
from datetime import datetime
import secrets
import string
import pytz
import os
import os.path
import re
import requests
import pathlib
from pathlib import Path
import pickle
# == Installed Modules ==
from alive_progress import alive_bar
import yaml
from importlib_resources import files
from Bio import SeqIO
from Bio.Data import IUPACData
import boto3
from botocore.exceptions import NoCredentialsError
from botocore import UNSIGNED
from botocore.config import Config

import pandas as pd
# == Project Modules ==


def compress_file(file_path: str):
	if not is_gzipped(file_path):
		# If not gzipped, compress the file
		with open(file_path, 'rb') as f_in, gzip.open(file_path + '.gz', 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
		print(f"File '{file_path}' compressed successfully.")
	if is_gzipped(file_path):
		cmd_rename = f"mv {file_path} {file_path}.gz"
		subprocess.run(cmd_rename, shell=True)
		print("This file is already compressed.")
		print(f"Created a copy of the file input on: {file_path}.gz")


def consolidate_s3_download(content, parent_folder):
	consolidated_downloadable = []
	for content_idx in range(0, len(content)):
		key = content[content_idx]['Key']
		if key != f"{parent_folder}/":
			consolidated_downloadable.append(content[content_idx])
	return consolidated_downloadable


def check_format(variable, data_type, paramater_name, default_value):
	if isinstance(variable, data_type):
		if variable == default_value:
			print(f"Please specify a value for the option --{paramater_name}.")
			exit(0)
		if data_type == str:
			if not check_iupac(variable):
				print(f"The string provided in --{paramater_name} is not a valid IUPAC representation: {variable}")
				exit(0)
		return data_type(variable)
	else:
		if data_type != str:
			try:
				data_type(variable)
				return data_type(variable)
			except ValueError:
				print(f"Invalid data type for {paramater_name}: '{variable}'. Please double-check the documentation.")
				exit(0)


def check_iupac(sequence):
	# Get the valid IUPAC characters (both unambiguous and ambiguous nucleotides)
	iupac_nucleotides = set(IUPACData.ambiguous_dna_letters)

	# Convert the sequence to uppercase and check if all characters are valid
	return all(base in iupac_nucleotides for base in sequence.upper())


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


def download_s3_objects(s3_bucket_name: str, s3_object_name: str, destination_path: str):
	"""
	Downloads a file or every file in a given AWS S3 bucket
	:param s3_bucket_name:
	:param s3_object_name:
	:param destination_path:
	:return:
	"""
	s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
	try:
		response = s3.list_objects_v2(Bucket=s3_bucket_name, Prefix=s3_object_name)
		content = response.get('Contents', [])
		clean_content = consolidate_s3_download(content, s3_object_name)  # Skip the S3 object that denotes the parent folder
	except NoCredentialsError:
		prRed("AWS issued credentials error. This is not an expected behavior. Please notify log this error on GitHub")
		exit(0)
	with alive_bar(len(clean_content), title=f'Downloading s3://{s3_bucket_name}/{s3_object_name}') as bar:
		for content_idx in range(0, len(clean_content)):
			key = clean_content[content_idx]['Key']
			source_file = key.split("/")[-1]
			destination_file = os.path.join(destination_path, source_file)
			if file_exists(destination_file):
				print(f"Skipping existing file: {destination_file}")
				bar()
				continue
			# == This downloads the data from the S3 Bucket without requiring AWS credentials
			pathlib.Path(destination_path).mkdir(parents=True, exist_ok=True)
			s3.download_file(s3_bucket_name, key, destination_file)
			bar()


def export_guides_by_editor(guide_df_by_editor_dict: dict, output_dir: (str, Path)):
	editors_list = []
	for editor in guide_df_by_editor_dict:
		editors_list.append(editor)
		guide_df = pd.DataFrame(guide_df_by_editor_dict[editor][0])
		filepath = f"{output_dir}{editor}.pkl"
		# Create output directory if non-existent
		set_export(output_dir)
		with open(filepath, 'wb') as guide_df_handle:
			pickle.dump(guide_df, guide_df_handle)
	return editors_list


def file_exists(file_path):
	return os.path.exists(file_path)


def group_guide_table(guide_table_path: pathlib.Path, editor_filter: (list, str)):
	guides_df = pd.read_csv(guide_table_path)
	if editor_filter:
		guides_df = guides_df[guides_df['Editor'].isin(editor_filter)]
	try:
		grouped_guides_df = guides_df.groupby('Editor')
	except KeyError:
		print("Column name 'Editor' not found in <Guides_found.csv>. Please check the file path and try again.")
		exit(0)
	editor_expanded_dictionary = {}
	for editor, stats in grouped_guides_df:
		editor_expanded_dictionary.setdefault(editor, []).append(stats)
	return editor_expanded_dictionary


def handle_shell_exception(subprocess_result, shell_command, verbose: bool):
	# === Silverplate Errors
	#   == File not found
	if re.search("FileNotFound", str(subprocess_result)):
		prRed(f"QUERY FILEPATH NOT FOUND: Please note:"
			  f"- This is the invalid filepath: {shell_command}"
			  f"- If providing more than one query path make sure the paths are comma-separated")
		exit(0)
	#	== Config Not Found
	if re.search("ConfigNotFound", str(subprocess_result)):
		prRed(f" Couldn't find the file {shell_command}."
			  f" Please double-check the path to <medit_database> and provide the path via '-d' option")
		exit(0)
	# === Handle SMK exceptions through subprocess
	#   == Unlock directory if necessary for SMK run
	if re.findall("Directory cannot be locked.", subprocess_result.stdout) or re.findall("Directory cannot be locked.", subprocess_result.stderr):
		print("--> Target directory locked. Unlocking...")
		unlock_smk_command = f"{shell_command} --unlock"
		launch_shell_cmd(unlock_smk_command, verbose)
		launch_shell_cmd(shell_command, verbose)
		return
	#	== Handle Missing Output Error in Snakemake
	if re.findall("MissingOutputException", subprocess_result.stdout):
		prRed("--> The association between Inputs/Outputs is not valid. Outputs may not be generated")
		return
	#   == Skipping rule call that has already been completed
	if re.findall(r"ValueError: min\(\) arg is an empty sequence", subprocess_result.stderr):
		print("--> A consensus FASTA has already been generated for this job. Skipping.")
		return
	if not re.findall(r"ValueError: min\(\) arg is an empty sequence", subprocess_result.stderr):
		if verbose:
			print(subprocess_result.stderr)
			prGreen(subprocess_result.stdout)
		return


def is_bgzipped(file_path: str) -> bool:
	with open(file_path, 'rb') as f:
		# Check if the file starts with the bgzip magic bytes
		# BGZF magic bytes: 1f 8b 08 04
		return f.read(4) == b'\x1f\x8b\x08\x04'


def is_gzipped(file_path: str):
	with open(file_path, 'rb') as f:
		# Check if the file starts with the gzip magic bytes
		return f.read(2) == b'\x1f\x8b'


def launch_shell_cmd(command: str, verbose=False, **kwargs):
	message = kwargs.get('message', False)
	check_exist = kwargs.get('check_exist', False)
	if message:
		verbose = False
		print(message)
	if check_exist:
		if os.path.isfile(check_exist):
			print(f"File {check_exist} exists. Skipping process.")
			return
	if verbose:
		prCyan(f"--> Invoking command-line call:\n{command}")

	result = subprocess.run(command,
	                        shell=True,
	                        stderr=subprocess.PIPE,
	                        stdout=subprocess.PIPE,
							text=True,
	                        universal_newlines=True
	                        )
	handle_shell_exception(result, command, verbose)


def list_files_by_extension(root_path, extension: str):
	file_list = []
	for root, dirs, files in os.walk(root_path, topdown=False):
		for name in files:
			if name.endswith(extension):
				file_list.append(os.path.join(root, name))
	return file_list


def offtarget_mode_formatting(mode, dynamic_config_guidepred):
	"""
	This function retrieves the relevant sequence IDs of VCFs and reference
	genome to help build the offtarget configuration file
	:param mode: current mEdit run mode
	:param dynamic_config_guidepred: YAML-imported Config dict
	:return: Tuple of size n>=1 with the sequence IDs and respective labels: Either 'main_ref' or 'extended'
	"""
	sequence_id = [(dynamic_config_guidepred['sequence_id'][0]), 'main_ref']
	if mode != 'fast':
		sequence_id.extend((dynamic_config_guidepred['vcf_id'], 'extended_ref'))
	return sequence_id


def pickle_chromosomes(genome_fasta, output_dir):
	records = SeqIO.parse(open(genome_fasta, 'rt'), "fasta")
	with alive_bar(25, title=f'Serializing human chromosomes') as bar:
		for record in records:
			if re.search(r"chr\w{0,2}$", record.id, re.IGNORECASE):
				outfile = f"{output_dir}/{record.id}.pkl"
				with open(outfile, 'ab') as gfile:
					pickle.dump(record, gfile)
					bar()


def prCyan(skk):
	print("\033[96m {}\033[00m".format(skk))


def prGreen(skk):
	print("\033[92m {}\033[00m".format(skk))


def prRed(skk):
	print("\033[0;31;47m {}\033[00m".format(skk))


def project_file_path(path_from_toplevel: str, filename: str):
	"""
	There are two top-level directories in the current version of mEdit: snakemake and config
	From either of these paths, the respective *.smk and *.yaml files can be accessed
	:param path_from_toplevel:
	:param filename:
	:return:
	"""
	return str(files(path_from_toplevel).joinpath(filename))


def set_export(outdir: str):
	if os.path.exists(outdir):
		pass
		# print(f'--> Skipping directory creation: {outdir}')
	# Create outdir only if it doesn't exist
	if not os.path.exists(outdir):
		print(f'Directory created on: {outdir}')
		pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
	return outdir


def validate_editor_list(editor_request: str, built_in_editors: list, parameter_string: str):
	validated_editors = []
	editor_list = editor_request.split(",")
	for editor in editor_list:
		green_light = False
		for built_in_editor in built_in_editors:
			clean_editor= editor.strip()
			if re.search(clean_editor, built_in_editor, re.IGNORECASE):
				green_light = True
				validated_editors.append(built_in_editor)
		if not green_light:
			print(
				f"NEW RE_STRATEGY: At least one of the requested editors is currently not built into mEdit: '{editor}'.\n "
				f"Full Set of Refs: {built_in_editors}"
				f"Please call 'medit list --help' to see directions on how to obtain the current list of editors.\n"
				f"Alternatively, use the '{parameter_string} custom' parameter to customize your own editor for this run.\n")
			exit(0)
	return validated_editors


def write_yaml_to_file(py_obj, filename: str):
	with open(f'{filename}', 'w', ) as f:
		yaml.safe_dump(py_obj, f, sort_keys=False, default_style='"')
	print(f'--> Configuration file created: {filename}')

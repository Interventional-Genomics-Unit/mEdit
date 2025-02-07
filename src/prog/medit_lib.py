# == Native Modules ==
from datetime import datetime
import gzip
import itertools
import logging
import os
import os.path
import pathlib
from pathlib import Path
import pickle
import pytz
import re
import requests
import secrets
import select
import shutil
import string
import subprocess
from typing import NamedTuple
# == Installed Modules ==
from alive_progress import alive_bar
from importlib_resources import files
from Bio import SeqIO
from Bio.Data import IUPACData
import boto3
from botocore.exceptions import NoCredentialsError
from botocore import UNSIGNED
from botocore.config import Config
import pandas as pd
import yaml
# == Project Modules ==

# === Define Classes ===


class QuotedStringDumper(yaml.Dumper):
	pass


# Define a structured result type
class SubprocessResult(NamedTuple):
	stdout: str
	stderr: str
	exit_code: int


# === Define Functions


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


def combine_guide_tables(combined_guide_report,path_to_report,genome_type):
	# if guide report exists add to compiled list of guides to run in guidescan2
	if file_exists(path_to_report):
		df = pd.read_csv(path_to_report)
		if df.empty:
			pass
		elif genome_type == 'main_ref':
			df.insert(4, "Alt Guide Impact", [None] * df.shape[0])
			df.insert(5, "Alt Genome", [None] * df.shape[0])
			combined_guide_report = pd.concat([combined_guide_report,df])
		elif genome_type == 'extended':
			select_columns = [x for x in df.columns if "Ref " not in x]
			df = df.loc[df['Alt Guide Impact'] != "PAM changed & removed",select_columns]
			df.columns = df.columns.to_list()[:5] +[x.replace("Alt ","") for x in df.columns[5:]]
			combined_guide_report = pd.concat([combined_guide_report, df])
		else:
			pass
	return combined_guide_report


def dataframe_duplet_to_dict(df, colA, colB):
	"""
	Converts two columns of a pandas DataFrame into a dictionary.

	Parameters:
		df (pd.DataFrame): The input DataFrame.
		colA (str): The column to use as keys.
		colB (str): The column to use as values.

	Returns:
		dict: A dictionary with colA as keys and colB as values.
	"""
	if colA not in df.columns or colB not in df.columns:
		raise ValueError(f"Columns '{colA}' or '{colB}' are not in the DataFrame.")

	return dict(zip(df[colA], df[colB]))


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


def get_public_gdrive_items(folder_id):
	"""
	Retrieves a list of file IDs, names, and types (file/folder) from a public Google Drive folder.
	:param folder_id: Google Drive public folder ID
	:return: List of (file_id, file_name, is_folder)
	"""
	url = f"https://drive.google.com/drive/folders/{folder_id}"
	response = requests.get(url)

	if response.status_code != 200:
		print(f"Failed to access the folder {folder_id}. Ensure it's public.")
		return []

	# Extract file IDs, names, and check if they are folders (based on class name patterns)
	matches = re.findall(r'data-id="(.*?)".*?data-tooltip="(.*?)".*?class="(.*?)"', response.text)

	items = []
	for file_id, file_name, class_info in matches:
		is_folder = "folder" in class_info.lower()  # If "folder" appears in class, it's a folder
		items.append((file_id, file_name, is_folder))

	return items


def download_public_gdrive_file(file_id, file_name, destination_path):
	"""
	Downloads a publicly shared Google Drive file. If manual confirmation is needed,
	prompts the user to visit the link and provide the direct download URL.
	"""
	session = requests.Session()
	download_url = f"https://drive.google.com/uc?export=download&id={file_id}"

	response = session.get(download_url, stream=True)

	# If the file is too large, Google blocks direct download
	if "Google Drive can't scan this file for viruses." in response.text:
		print("\n🚨 Google Drive requires confirmation for large file downloads.")
		print("👉 Please open the following URL in your browser:")
		print(f"\n🔗 {download_url}\n")
		print("1️⃣ Click 'Download anyway' when prompted.")
		print("2️⃣ Copy the final download link (the one that actually starts downloading).")
		direct_url = input("\n🔽 Paste the final download link here: ").strip()
		if not direct_url:
			print("❌ No URL provided. Aborting.")
			return
		response = session.get(direct_url, stream=True)

	# Ensure target directory exists
	pathlib.Path(destination_path).mkdir(parents=True, exist_ok=True)
	file_path = os.path.join(destination_path, file_name)

	# Save the file
	with open(file_path, "wb") as f:
		for chunk in response.iter_content(1024 * 1024):  # Download in 1MB chunks
			if chunk:
				f.write(chunk)

	print(f"\n✅ Downloaded: {file_path}")


def download_gdrive_folder(folder_id, destination_path):
	"""
	Recursively downloads all files and subfolders from a public Google Drive folder.
	:param folder_id: Google Drive public folder ID
	:param destination_path: Local directory to save files
	"""
	items = get_public_gdrive_items(folder_id)

	if not items:
		print(f"No items found in folder: {folder_id}")
		return

	print(f"Found {len(items)} items in folder {folder_id}. Downloading...")

	with alive_bar(len(items), title=f"Downloading Google Drive Folder: {destination_path}") as bar:
		for file_id, file_name, is_folder in items:
			file_name = re.sub(r"^(Compressed archive: |Text: |Spreadsheet: |Presentation: )", "", file_name).strip()
			local_path = os.path.join(destination_path, file_name)

			if is_folder:
				print(f"\nEntering subfolder: {file_name}")
				download_gdrive_folder(file_id, local_path)  # Recursive call for subfolders
			else:
				download_public_gdrive_file(file_id, file_name, destination_path)

			bar()

	print(f"Finished downloading folder: {destination_path}")


def download_s3_objects(s3_bucket_name: str, s3_object_name: str, destination_path: str):
	"""
	Downloads a file or every file in a given AWS S3 bucket
	:param s3_bucket_name:
	:param s3_object_name:
	:param destination_path:
	:return: skip_flag: True if the download has been skipped, False otherwise
	"""
	skip_flag = False
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
				skip_flag = True
				continue
			# == This downloads the data from the S3 Bucket without requiring AWS credentials
			pathlib.Path(destination_path).mkdir(parents=True, exist_ok=True)
			s3.download_file(s3_bucket_name, key, destination_file)
			bar()
	return skip_flag


def expand_pam(pam, mode='ref'):
	# guidescan will not accept anything aside ACGTN
	# if a pam has a different base this will make a list of variants to be give to guidescan
	iupac_codes = {"A": "A", "T": "T", "C": "C",
				   "G": "G", "N": "N", "R": "AG", "Y": "CT",
				   "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
				   "D": "AGT", "H": "ACT", "V": "ACG"}
	pam = pam.upper()
	try:
		expanded = [list(iupac_codes[base]) for base in pam]
		expanded_pams = [''.join(x) for x in itertools.product(*expanded)]
		input_file_pam = expanded_pams[0]
		cmd_line_alt_pams = ",".join(expanded_pams[1:]) if len(expanded_pams) > 1 else "no_alt_pam"
	except KeyError:
		print(f"invalid PAM given: {pam}")
		print(f"PAM accepted characters 'ATCGNRYSWKMBDHV' ")
		exit(0)

	return [input_file_pam, cmd_line_alt_pams]


def export_guides_by_editor(guide_df_by_editor_dict: dict, output_dir: (str, Path)):
	editors_list = []
	set_export("/".join(str(output_dir).split("/")[:-1]))
	for editor in guide_df_by_editor_dict:
		editors_list.append(editor)
		guide_df = pd.DataFrame(guide_df_by_editor_dict[editor][0])
		filepath = f"{output_dir}{editor}.pkl"
		# Create output directory if non-existent
		with open(filepath, 'wb') as guide_df_handle:
			pickle.dump(guide_df, guide_df_handle)
	return editors_list


def file_exists(file_path):
	return os.path.exists(file_path)


def group_guide_table(allguides_df, editor_filter: (list, str)):
	if editor_filter:
		allguides_df = allguides_df[allguides_df['Editor'].isin(editor_filter)]
	try:
		grouped_guides_df = allguides_df.groupby('Editor')
	except KeyError:
		prRed("The column name format found in <Guides_found.csv> is not compliant with this version of mEdit."
			  " Please check the file path and try again.")
		exit(0)
	editor_expanded_dictionary = {}
	for editor, stats in grouped_guides_df:
		editor_expanded_dictionary.setdefault(editor, []).append(stats)
	return editor_expanded_dictionary


def handle_offtarget_request(mode, meditdb_path, requested_genomes):
	gscan_indices_path = f"{meditdb_path}/gscan_indices"
	issue_warning = False
	n_genomes_found = []
	genome_ids = [t[0] for t in requested_genomes]
	n_genomes_to_process = len(requested_genomes)

	if mode != 'fast':
		if is_empty_or_nonexistent(gscan_indices_path):
			issue_warning =True
		elif not is_empty_or_nonexistent(gscan_indices_path):
			# Precompute regex patterns for all requested genomes (case-insensitive, whole-word match)
			regex_patterns = [
				re.compile(rf"{genome}.*\.index\..*", re.IGNORECASE)  # Match whole words only
				for genome in genome_ids
			]

			for root, dirs, files in os.walk(gscan_indices_path):
				for file in files:
					# Check if the file matches ANY requested genome
					for pattern in regex_patterns:
						if pattern.search(file):
							# Extract the matched genome name (case-insensitive)
							matched_genome = genome_ids[regex_patterns.index(pattern)]
							n_genomes_found.append(matched_genome.lower())  # Case-insensitive uniqueness
							break  # Stop checking other genomes for this file

			if len(requested_genomes) > len(n_genomes_found):
				n_genomes_to_process = len(requested_genomes) - (len(n_genomes_found) + 1)
				issue_warning = True

	if issue_warning:
		prRed(f"!! WARNING !! Read this carefully:\n"
			  f"- This execution is about to carry out an off-target analysis "
			  f"utilizing more than one genome.\n <{genome_ids}>"
			  f"- The current mEdit database does not contain the necessary support files on this path:\n"
			  f"[{gscan_indices_path}]\n"
			  f"- You'll need {n_genomes_to_process} genomes processed for the current run setup.\n "
			  f"- mEdit can generate those files, but it is highly recommended to do that in a multi-processor\n"
			  f"environment (such as an HPC) that can be leveraged to distribute the processes.\n"
			  f"- Otherwise, this process can take up to several days to complete depending on how many "
			  f"genomes are processed in series.")

		# Prompt the user to type 'y' to continue
		user_input = input("Would you like to continue? (y/n): ").strip().lower()

		if user_input != 'y':
			user_input = input(
				"Would you like to switch to 'fast' mode and proceed with Hg38 genome only? (y/n): ").strip().lower()
			if user_input == 'y':
				return True, 'fast'
			print("Operation aborted by the user.")
			return False, mode

	# Continue with the off-target analysis
	print("Resuming off-target analysis...")
	return True, mode


def handle_shell_exception(result: SubprocessResult, shell_command: str, verbose: bool) -> None:
	"""Handle errors in subprocess output, including Snakemake locks."""
	# === Silverplate Errors
	#   == File not found
	if re.search("FileNotFound", result.stderr):
		prRed(f"QUERY FILEPATH NOT FOUND: Invalid filepath in command: {shell_command}")
		exit(1)

	# === Snakemake-Specific Errors
	#   == Unlock directory if locked
	if "Directory cannot be locked." in result.stdout + result.stderr:
		prCyan("--> Target directory locked. Unlocking...")
		unlock_command = f"{shell_command} --unlock"  # Replace with actual unlock command
		unlock_result = launch_shell_cmd(unlock_command, verbose)

		if unlock_result.exit_code != 0:
			prRed("--> Failed to unlock directory!")
			exit(1)

		# Retry original command after unlocking
		prCyan("--> Retrying original command...")
		retry_result = launch_shell_cmd(shell_command, verbose)
		handle_shell_exception(retry_result, shell_command, verbose)
		return

	#   == Other Snakemake errors
	if "MissingOutputException" in result.stdout:
		prRed("--> Missing outputs. Check input/output associations.")
		exit(1)

	# === General errors
	if result.exit_code != 0:
		prRed(f"Command failed with exit code {result.exit_code}: {shell_command}")
		exit(1)


def init_logging(logfile_path):
	# Configure the logging system
	logging.basicConfig(
		level=logging.DEBUG,  # Set the minimum log level (DEBUG logs everything)
		format='%(asctime)s [%(levelname)s] %(message)s',  # Define log format
		handlers=[
			logging.FileHandler(logfile_path),  # Log to file
		]
	)


def is_bgzipped(file_path: str) -> bool:
	with open(file_path, 'rb') as f:
		# Check if the file starts with the bgzip magic bytes
		# BGZF magic bytes: 1f 8b 08 04
		return f.read(4) == b'\x1f\x8b\x08\x04'


def is_gzipped(file_path: str):
	with open(file_path, 'rb') as f:
		# Check if the file starts with the gzip magic bytes
		return f.read(2) == b'\x1f\x8b'


# def launch_shell_cmd(command: str, verbose=False, **kwargs):
# 	message = kwargs.get('message', False)
# 	check_exist = kwargs.get('check_exist', False)
#
# 	command = f"stdbuf -oL -eL {command}"
#
# 	if message:
# 		verbose = False
# 		print(message)
# 	if check_exist:
# 		if os.path.isfile(check_exist):
# 			print(f"File {check_exist} exists. Skipping process.")
# 			return
# 	if verbose:
# 		prCyan(f"--> Invoking command-line call:\n{command}")
#
# 	# CURRENTLY IN REVISION
# 	process = subprocess.Popen(
# 		command,
# 		shell=True,
# 		stdout=subprocess.PIPE,
# 		stderr=subprocess.PIPE,
# 		text=True,  # Ensures strings are returned instead of bytes
# 		universal_newlines=True
# 	)
#
# 	# Create a polling object to monitor both streams
# 	poll = select.poll()
# 	poll.register(process.stdout.fileno(), select.POLLIN)
# 	poll.register(process.stderr.fileno(), select.POLLIN)
#
# 	try:
# 		while True:
# 			# Poll both streams
# 			events = poll.poll()
# 			for fd, event in events:
# 				if fd == process.stdout.fileno():
# 					# Read and print stdout
# 					line = process.stdout.readline()
# 					if line:
# 						prGreen(line, False)
# 				elif fd == process.stderr.fileno():
# 					# Read and print stderr
# 					line = process.stderr.readline()
# 					if line:
# 						prCyan(line, False)
#
# 			# Exit the loop if the process finishes
# 			if process.poll() is not None:
# 				break
#
# 		# Ensure remaining lines are flushed
# 		for line in process.stdout:
# 			prRed(line, False)
# 		for line in process.stderr:
# 			prRed(line, False)
#
# 	except Exception as e:
# 		print(f"An error occurred: {e}")
# 		process.terminate()
# 		process.wait()


def is_empty_or_nonexistent(directory_path):
	"""
	Checks if a directory is either empty or does not exist.

	:param directory_path: Path to the directory.
	:return: True if the directory does not exist or is empty, False otherwise.
	"""
	return not os.path.isdir(directory_path) or not any(os.scandir(directory_path))


def launch_shell_cmd(command: str, verbose: bool = False, **kwargs) -> SubprocessResult:
	"""Execute a shell command with real-time output capture and error handling."""
	message = kwargs.get('message')
	check_exist = kwargs.get('check_exist')

	# Pre-command checks
	if check_exist and os.path.exists(check_exist):
		print(f"File {check_exist} exists. Skipping.")
		return SubprocessResult("", "", 0)

	# Command modifications
	if 'decompress' in command:
		command = command.replace('decompress', 'pigz -p $(nproc)' if shutil.which("pigz") else 'gzip')

	command = f"stdbuf -oL -eL {command}"  # Line-buffered output

	if message:
		print(message)

	if verbose:
		prCyan(f"--> Executing:\n{command}")

	# Run the command
	process = subprocess.Popen(
		command,
		shell=True,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,
		text=True,
		universal_newlines=True
	)

	# Capture output in real-time
	stdout_lines = []
	stderr_lines = []
	poller = select.poll()
	poller.register(process.stdout, select.POLLIN)
	poller.register(process.stderr, select.POLLIN)

	try:
		while True:
			events = poller.poll()
			for fd, event in events:
				if fd == process.stdout.fileno():
					line = process.stdout.readline()
					stdout_lines.append(line)
					if verbose:
						prGreen(line.strip())
				elif fd == process.stderr.fileno():
					line = process.stderr.readline()
					stderr_lines.append(line)
					if verbose:
						print(line.strip())

			if process.poll() is not None:
				break

		# Capture any remaining output
		stdout_remainder = process.stdout.read()
		stderr_remainder = process.stderr.read()
		stdout_lines.append(stdout_remainder)
		stderr_lines.append(stderr_remainder)

	except Exception as e:
		process.terminate()
		return SubprocessResult("".join(stdout_lines), "".join(stderr_lines), 1)

	return SubprocessResult(
		stdout="".join(stdout_lines),
		stderr="".join(stderr_lines),
		exit_code=process.returncode
	)


def list_files_by_extension(root_path, extension: str):
	file_list = []
	for root, dirs, files in os.walk(root_path, topdown=False):
		for name in files:
			if name.endswith(extension):
				file_list.append(os.path.join(root, name))
	return file_list


def offtarget_mode_formatting(mode, reference_genome, dynamic_config_guidepred):
	"""
	This function retrieves the relevant sequence IDs of VCFs and reference
	genome to help build the offtarget configuration file
	:param reference_genome: Identifier string of the reference genome
	:param mode: current mEdit run mode
	:param dynamic_config_guidepred: YAML-imported Config dict
	:return: Tuple of size n>=1 with the sequence IDs and respective labels: Either 'main_ref' or 'extended'
	"""
	sequence_id = [(reference_genome, 'main_ref')]
	if mode != 'fast':
		for vcf_id in dynamic_config_guidepred['vcf_id']:
			sequence_id.extend([(vcf_id, 'extended')])
	return sequence_id


def prCyan(skk, newline=True):
	if newline:
		print("\033[96m {}\033[00m".format(skk))
	else:
		print("\033[96m {}\033[00m".format(skk), end='')


def prGreen(skk, newline=True):
	if newline:
		print("\033[92m {}\033[00m".format(skk))
	else:
		print("\033[92m {}\033[00m".format(skk), end='')


def prRed(skk, newline=True):
	if newline:
		print("\033[0;31;47m {}\033[00m".format(skk))
	else:
		print("\033[0;31;47m {}\033[00m".format(skk), end='')


def prYellow(skk, newline=True):
	if newline:
		print("\033[93m {}\033[00m".format(skk))
	else:
		print("\033[93m {}\033[00m".format(skk), end='')


def project_file_path(path_from_toplevel: str, filename: str):
	"""
	There are two top-level directories in the current version of mEdit: snakemake and config
	From either of these paths, the respective *.smk and *.yaml files can be accessed
	:param path_from_toplevel:
	:param filename:
	:return:
	"""
	return str(files(path_from_toplevel).joinpath(filename))


def quoted_string_representer(dumper, data):
	return dumper.represent_scalar('tag:yaml.org,2002:str', data, style="'")


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
				f"Full Set of Refs: {built_in_editors}"
				f"Please call 'medit list --help' to see directions on how to obtain the current list of editors.\n"
				f"Alternatively, use the '{parameter_string} custom' parameter to customize your own editor for this run.\n")
			exit(0)
	return validated_editors


def write_yaml_to_file(py_obj, filename: str):
	# Add the custom representer to ensure strings are quoted
	yaml.add_representer(str, quoted_string_representer, Dumper=QuotedStringDumper)
	with open(f'{filename}', 'w', ) as f:
		yaml.safe_dump(py_obj, f, sort_keys=False, default_style='"')
		# yaml.dump(py_obj, f, sort_keys=False, Dumper=QuotedStringDumper, default_flow_style=False)
	print(f'--> Configuration file created: {filename}')

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
import pathlib
import pickle
# == Installed Modules ==
from alive_progress import alive_bar
import yaml
from importlib_resources import files
from Bio import SeqIO
import boto3
from botocore.exceptions import NoCredentialsError
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
	s3 = boto3.client('s3')
	try:
		response = s3.list_objects_v2(Bucket=s3_bucket_name, Prefix=s3_object_name)
		content = response.get('Contents', [])
		clean_content = consolidate_s3_download(content, s3_object_name)  # Skip the S3 object that denotes the parent folder
	except NoCredentialsError:
		prRed("AWS Credentials not available. Please sign in then try again!")
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
			pathlib.Path(destination_path).mkdir(parents=True, exist_ok=True)
			s3.download_file(s3_bucket_name, key, destination_file)
			bar()


def file_exists(file_path):
	return os.path.exists(file_path)


def handle_shell_exception(subprocess_result, shell_command, verbose: bool):
	# === Handle SMK exceptions through subprocess
	#   == Unlock directory if necessary for SMK run
	if re.findall("Directory cannot be locked.", subprocess_result.stdout):
		print("--> Target directory locked. Unlocking...")
		unlock_smk_command = f"{shell_command} --unlock"
		launch_shell_cmd(unlock_smk_command, verbose)
		launch_shell_cmd(shell_command, verbose)
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


def parse_editor_request(request):
	processed_request = request
	if type(processed_request) is list:
		processed_request = request.split(',')

	return processed_request


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
		print(f'--> Skipping directory creation: {outdir}')
	# Create outdir only if it doesn't exist
	if not os.path.exists(outdir):
		print(f'Directory created on: {outdir}')
		pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
	return outdir


def write_yaml_to_file(py_obj, filename: str):
	with open(f'{filename}', 'w', ) as f:
		yaml.safe_dump(py_obj, f, sort_keys=False, default_style='"')
	print(f'--> Configuration file created: {filename}')

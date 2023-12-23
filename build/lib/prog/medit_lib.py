# == Native Modules ==
import subprocess
import gzip
import shutil
from datetime import datetime
import secrets
import string
import pytz
import os
import re
import pickle
# == Installed Modules ==
import yaml
from Bio import SeqIO
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
		print("This VCF file is already compressed.")
		print(f"Created a copy of the VCF file input on: {file_path}.gz")


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


def is_gzipped(file_path: str):
	with open(file_path, 'rb') as f:
		# Check if the file starts with the gzip magic bytes
		return f.read(2) == b'\x1f\x8b'


def launch_shell_cmd(command: str, verbose: bool):
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
	for record in records:
		if re.search(r"chr\w{0,2}$", record.id):
			outfile = f"{output_dir}/{record.id}.pkl"
			with open(outfile, 'ab') as gfile:
				print(f"Serializing chromosome {record.id}")
				pickle.dump(record, gfile)


def prCyan(skk):
	print("\033[96m {}\033[00m" .format(skk))


def prGreen(skk):
	print("\033[92m {}\033[00m" .format(skk))


def set_export(outdir: str):
	if os.path.exists(outdir):
		print(f'--> Skipping directory creation: {outdir}')
	# Create outdir only if it doesn't exist
	if not os.path.exists(outdir):
		print(f'Directory created on: {outdir}')
		os.makedirs(outdir)
	return outdir


def write_yaml_to_file(py_obj, filename: str):
	with open(f'{filename}', 'w',) as f:
		yaml.safe_dump(py_obj, f, sort_keys=False, default_style='"')
	print(f'--> Configuration file created: {filename}')

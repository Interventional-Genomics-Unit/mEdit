# == Native Modules ==
import subprocess
import gzip
import shutil
from datetime import datetime
import secrets
import string
import pytz
import os
# == Installed Modules ==
import yaml
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


def is_gzipped(file_path: str):
	with open(file_path, 'rb') as f:
		# Check if the file starts with the gzip magic bytes
		return f.read(2) == b'\x1f\x8b'


def launch_shell_cmd(command: str):
	print(f"Invoking command-line call:\n{command}")
	subprocess.run(command, shell=True)


def set_export(outdir: str):
	# Create outdir and all missing parent directories
	os.makedirs(outdir, exist_ok=True)
	return outdir


def write_yaml_to_file(py_obj, filename: str):
	with open(f'{filename}', 'w',) as f:
		yaml.safe_dump(py_obj, f, sort_keys=False, default_style='"')
	print(f'Configuration file sucessfully written to: {filename}')

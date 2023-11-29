# == Native Modules ==
from datetime import datetime
import secrets
import string
import pytz
# == Installed Modules ==
import yaml
# == Project Modules ==
from programs.guide_prediction import run as guide_prediction
from programs.db_set import run as db_set
from programs.arguments import run as parse_arguments


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


def write_yaml_to_file(py_obj, filename):
	with open(f'{filename}.yaml', 'w',) as f:
		yaml.dump(py_obj, f, sort_keys=False)
	print('Written to file successfully')


def main():
	# === Call argument parsing function ===
	args = parse_arguments()
	# mEdit Program
	program = args.program

	# Assign jobtag and run mode to config
	jobtag = date_tag()

	if program == "guide_prediction":
		guide_prediction(args, jobtag)

	# == Database Parameters
	if program == "db_set":
		db_set(args, jobtag)


if __name__ == "__main__":
	main()

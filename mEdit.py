# == Native Modules ==
from datetime import datetime
import secrets
import string
import pytz
# == Installed Modules ==
import yaml
# == Project Modules ==
from programs.guide_prediction import guide_prediction as guide_prediction
from programs.db_set import dbset as db_set
from programs.arguments import parse_arguments as parse_arguments
from programs.medit_lib import date_tag as date_tag


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

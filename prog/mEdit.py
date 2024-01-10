# == Native Modules ==
import os
# == Installed Modules ==
# == Project Modules ==
from prog.guide_prediction import guide_prediction as guide_prediction
from prog.db_set import dbset as db_set
from prog.arguments import parse_arguments
from prog.medit_lib import date_tag


def main():
	# === Call argument parsing function ===
	args = parse_arguments()
	# mEdit Program
	program = args.program
	# Assign jobtag and run mode to config
	args.user_jobtag = True
	try:
		jobtag = args.jobtag
		if not jobtag:
			jobtag = date_tag()
	except AttributeError:
		jobtag = date_tag()
		args.user_jobtag = False
	# Run mEdit Programs
	if program == "guide_prediction":
		guide_prediction(args, jobtag)

	# == Database Parameters
	if program == "db_set":
		db_set(args)


if __name__ == "__main__":
	main()

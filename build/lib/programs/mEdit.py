# == Native Modules ==
# == Installed Modules ==
# == Project Modules ==
from prog.guide_prediction import guide_prediction as guide_prediction
from db_set import dbset as db_set
from arguments import parse_arguments
from medit_lib import date_tag


def main():
	# === Call argument parsing function ===
	args = parse_arguments()
	# mEdit Program
	program = args.program
	jobtag = args.jobtag
	args.user_jobtag = True
	# Assign jobtag and run mode to config
	if not jobtag:
		jobtag = date_tag()
		args.user_jobtag = False

	if program == "guide_prediction":
		guide_prediction(args, jobtag)

	# == Database Parameters
	if program == "db_set":
		db_set(args)


if __name__ == "__main__":
	main()

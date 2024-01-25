# == Native Modules
import pickle
from os.path import abspath
# == Installed Modules
import yaml
# == Project Modules


def ls_editors(args):
	# == Load Run Parameters values ==
	print_editors = args.editors
	print_base_editors = args.base_editors

	# === Load Database Path ===
	db_path_full = f"{abspath(args.db_path)}/medit_database"
	config_db_path = f"{db_path_full}/config_db/config_db.yaml"

	# === Load configuration file ===
	with open(config_db_path, 'r') as config_handle:
		config_db_obj = yaml.safe_load(config_handle)

	# === Load Editors Lists From Path Specified on config_db.yaml ===
	if print_editors:
		with open(str(config_db_obj["editors"]), 'rb') as editors_pkl_handle:
			editors_dict = pickle.load(editors_pkl_handle)
			print(f"Available editors:\n {list(editors_dict['clinical'].keys())}")
	if print_base_editors:
		with open(str(config_db_obj["base_editors"]), 'rb') as be_pkl_handle:
			base_editors_dict = pickle.load(be_pkl_handle)
			print(f"Available base editors:\n {list(base_editors_dict['all'].keys())}")

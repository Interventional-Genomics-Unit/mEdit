# == Native Modules
import pickle
from os.path import abspath, isdir
# == Installed Modules
import yaml
# == Project Modules
from prog.medit_lib import file_exists


def ls_editors(args):
	# == Load Run Parameters values ==
	# print_editors = args.editors
	# print_base_editors = args.base_editors

	# === Load Database Path ===
	db_path_full = f"{abspath(args.db_path)}/medit_database"
	config_db_path = f"{db_path_full}/config_db/config_db.yaml"

	if not file_exists(db_path_full):
		print("The database path directory could not be found.")
		exit(0)

	# === Load configuration file ===
	with open(config_db_path, 'r') as config_handle:
		config_db_obj = yaml.safe_load(config_handle)

	# === Load Editors Lists From Path Specified on config_db.yaml ===
	# if print_editors:
	with open(str(config_db_obj["editors"]), 'rb') as editors_pkl_handle:
		editors_dict = pickle.load(editors_pkl_handle)
		print(f"Available endonuclease editors: ")
		for k,v in editors_dict['all'].items():
			print("-----------------------------")
			print(f"name: {k}")
			print(f"pam, pam_is_first: {v[0]}, {v[1]}")
			print(f"guide_len: {str(v[2])}")
			if v[1] == True:# show both staggered positions for 5' pam editors
				print(f"dsb_position: {v[3]},{v[3]+5}")
			else:
				print(f"dsb_position: {v[3]}")
			if len(v) > 4:
				print(f"notes: {v[4]}")
			print(f"5'-{v[0] if v[1] else ''}{'x' * int(v[2])}{'' if v[1] else v[0]}-3'")


	# if print_base_editors:
	with open(str(config_db_obj["base_editors"]), 'rb') as be_pkl_handle:
		base_editors_dict = pickle.load(be_pkl_handle)
		print(f"Available base editors: ")
		for k,v in base_editors_dict['all'].items():

			print("-----------------------------")
			print(f"name: {k}")
			print(f"pam, pam_is_first: {v[0][0]}, {v[0][1]}")
			print(f"guide_len: {str(v[0][2])}")
			print(f"edit_win: {str(v[0][3])[1:-1]}")
			print(f'target_base, result_base: {v[1][0][0]} ---> {v[1][0][1]}')
			if len(v[0]) > 4:
				print(f"notes: {v[0][4]}")
			print(f"5'-{v[0][0] if v[0][1] else ''}{'x' * int(v[0][2])}{'' if v[0][1] else v[0][0]}-3'")


# == Native Modules
from os.path import abspath
import pickle
# == Installed Modules
import yaml
# == Project Modules
from prog.medit_lib import group_guide_table


def offtarget_prediction(args, jobtag):
	# == Load Run Parameters values ==
	user_jobtag = args.user_jobtag
	root_dir = abspath(args.output)
	# == Set export paths tied to the SMK pipeline ==
	config_dir_path = f"{root_dir}/config"
	# == Set export paths for dynamic YAML files ==
	dynamic_config_path = f"{config_dir_path}/config_{jobtag}.yaml"

	# ->=== CONFIG FILES IMPORT ===<-
	with open(dynamic_config_path, 'r') as config_handle:
		config_template = yaml.safe_load(config_handle)
	# TODO: Create a config template for offtargets and Pull here to assign variables

	# === Import Variables from Configuration File ===
	run_name = config_template['run_name']
	mode = config_template['processing_mode']
	root_dir = config_template['output_directory']
	sequence_id = config_template['sequence_id']

	# == Set output paths ==
	guides_per_editor_path = f"{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/guides_per_editor.pkl"

	# === Recover Guide Prediction filepath
	guides_report_path = (f"{root_dir}/{mode}/jobs/{run_name}/"
						 f"guide_prediction-{sequence_id}/guides_report_ref/Guides_found.csv"),
	guide_search_params = (f"{root_dir}/{mode}/jobs/{run_name}/"
						   f"guide_prediction-{sequence_id}/dynamic_params/guide_search_params.pkl")

	grouped_guide_dict = group_guide_table(guides_report_path)
	with open(guides_per_editor_path, 'ab') as guides_per_editor_handle:
		pickle.dump(grouped_guide_dict, guides_per_editor_handle)


if __name__ == "__main__":
	main()

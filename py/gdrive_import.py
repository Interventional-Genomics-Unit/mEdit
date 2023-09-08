# Native modules
import os
import re
# Installed modules
import yaml
import gdown
import pandas as pd


def get_col_as_list(csv_path, target_col):
	list_of_links = pd.read_csv(csv_path)[target_col].tolist()
	return list_of_links


def get_excel_from_google_drive(url_list, id_list, start_tab_index=1):
	association_dict = {}

	for idx in range(len(url_list)):
		sheets_dict = {}
		url = url_list[idx]
		output_dir = get_absolute_path("jobs")
		output_file = f'{output_dir}temp.xlsx'

		gdown.download(url, output_file, quiet=False, fuzzy=True)  # Download the file using gdown
		xl = pd.read_excel(output_file, sheet_name=None, engine='openpyxl')
		sheet_names = list(xl.keys())[start_tab_index:]

		for sheet_name in sheet_names:
			cut_sheet_name = re.sub('_table', '', sheet_name, flags=re.IGNORECASE)
			sheet_split = re.split('_', cut_sheet_name)
			sheet_name_loop = ''
			for field in range(0, 3):
				try:
					sheet_name_loop += f"{sheet_split[field]}_"
				except IndexError:
					break
			clean_sheet_name = sheet_name_loop.strip('_')

			# Load the Excel file into a pandas DataFrame
			df = pd.read_excel(output_file, sheet_name=sheet_name)
			sheets_dict.setdefault(clean_sheet_name, df)

		# Remove Temps
		os.remove(output_file)

		clean_id = re.sub(r'(\S+) \S+', r'\1', id_list[idx])
		print(clean_id)
		association_dict[clean_id] = sheets_dict
	# Return the dictionary
	return association_dict


# Load config_render file
with open("db_interaction.yaml", "r") as f:
	config = yaml.safe_load(f)


def main():


	# Consolidate Cloud links and identifiers
	wiki_links_list = get_col_as_list(config["wiki_manifest_path"], config["links_column"])
	uniq_id_list = get_col_as_list(config["wiki_manifest_path"], config["unique_id_col"])

	# Download forms from Google-drive
	id_to_sheets_dict = get_excel_from_google_drive(wiki_links_list, uniq_id_list)


if __name__ == "__main__":
	main()

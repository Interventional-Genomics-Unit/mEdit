# == Native Modules
import pickle


def create_guidescan_infile(casoff_input_path, guides_seq, pam, guide_names, coords):
	# create input file for guidescan
	header = "id,sequence,pam,chromosome,position,sense"

	with open(casoff_input_path, 'w') as f:
		f.writelines(header + "\n")

		for guide_name, guide_seq, coord in zip(guide_names, guides_seq, coords):
			f.write(",".join([guide_name,
							  guide_seq.upper(),
							  pam.upper(),
							  f"chr{coord.split(':')[0]}",
							  coord.split(':')[1][:-1].split("-")[0],
							  coord[-1]+"\n"]))


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	guides_report_per_editor_path = str(snakemake.input.guides_per_editor_path)
	# === Outputs ===
	casoff_input_path = str(snakemake.output.casoff_input)
	# === Params ===
	pam_per_editor = str(snakemake.params.pam_per_editor_dict)

	# === Import Per-editor data sets
	with open(guides_report_per_editor_path, 'rb') as f:
		guides_report_per_editor = pickle.load(f)

	# Gather data from the imported pickle
	guides_seq = list(guides_report_per_editor.gRNA)
	guide_names = list(guides_report_per_editor.Guide_ID)
	coords = list(guides_report_per_editor.Coordinates + guides_report_per_editor.Strand)

	create_guidescan_infile(casoff_input_path, guides_seq, pam_per_editor, guide_names, coords)


if __name__ == "__main__":
	main()



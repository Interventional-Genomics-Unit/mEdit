# == Native Modules
import pickle

def create_guidescan_infile(casoff_input_path, guides,editor_pam, gnames,coords):
	## create input file for guidescan
	header = "id,sequence,pam,chromosome,position,sense"

	with open(casoff_input_path, 'w') as f:
		f.writelines(header + "\n")

		for gname,grna,coord  in zip(gnames,guides,coords):
			f.write(",".join([gname,
					grna.upper(),
					editor_pam.upper(),
					f"chr{coord.split(':')[0]}",
					coord.split(':')[1][:-1],coord[-1],"\n"]))

def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	guides_report_per_editor_path = str(snakemake.input.guides_per_editor_path)

	# === Outputs ===
	casoff_input_path = str(snakemake.output.casoff_input)

	# === Wildcards ===
	editor_pam = str(snakemake.input.editor_pam)

	# === Guide search params ===
	guides_report_per_editor = pickle.load(open(guides_report_per_editor_path,'rb'))

	guides,gnames, coords = list(guides_report_per_editor.gRNA),list(guides_report_per_editor.Guide_ID), list(guides_report_per_editor.Coordinates + guides_report_per_editor.Strand)

	create_guidescan_infile(casoff_input_path, guides, editor_pam, gnames, coords)


if __name__ == "__main__":
	main()
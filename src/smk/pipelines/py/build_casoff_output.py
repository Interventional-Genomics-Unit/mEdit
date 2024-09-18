# == Native Modules
import pickle

# == Installed Modules
import pandas as pd


# == Project Modules


def write_out_res(concat_offtarget_scores_dict,
				  concat_guides_report,
				  casoff_params,
				  offtarget_scores_out,
				  formatted_casoff_out):
	'''
	Adds two new columns to guides found table
	Creates and Aggregate summary text output
	'''
	df = pd.DataFrame(concat_offtarget_scores_dict)  # .rename(columns = {'0':'Off Target Score'})
	# this is just a precautionary if this script was run more than once
	if 'Num Off Targets per MM' in concat_guides_report.columns:
		concat_guides_report = concat_guides_report.drop(columns='Num Off Targets per MM')
	if 'Off Target Score' in concat_guides_report.columns:
		concat_guides_report = concat_guides_report.drop(columns='Off Target Score')

	# Update main guides table with specifity Scoress
	spec_scores = df.T.sort_index()['spec_score']
	df = df.drop(index='spec_score')
	df = df.astype('int')
	# Add off_targets summary
	df.index = pd.MultiIndex.from_tuples(list(df.index), names=['BulgeType', 'Number of Mismatches'])
	df = df.reset_index()
	ot_per_mm = df.loc[df['BulgeType'] == 'X']
	ot_per_mm = ot_per_mm.iloc[:, 2:].T.sort_index()
	ot_per_mm = ot_per_mm.stack().astype('str').groupby(level=0).apply('|'.join)
	ot_per_mm = ot_per_mm.rename('Guide_ID')  # ('Num Off Targets per Mismatch')
	concat_guides_report = concat_guides_report.sort_values('Guide_ID')
	concat_guides_report['Off Target Score'] = list(spec_scores)
	concat_guides_report['Num Off Targets per MM'] = list(ot_per_mm)

	# create off_target summary of totals
	if casoff_params[1] == 0:
		df = df.loc[df.BulgeType != 'RNA']
	if casoff_params[2] == 0:
		df = df.loc[df.BulgeType != 'DNA']
	offtarget_summary_report = df.pivot_table(columns=['BulgeType', 'Number of Mismatches'], aggfunc="sum")

	# writeout
	# sum_out = resultsfolder + 'OffTarget_Summary.txt'
	offtarget_summary_report.to_csv(offtarget_scores_out, sep='\t')

	concat_guides_report.to_csv(formatted_casoff_out, index=False)


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	offtarget_scores_list = list(snakemake.output.offtarget_scores)
	guides_report_per_editor_list = list(snakemake.input.guides_per_editor_path)
	# === Outputs ===
	formatted_casoff_out = str(snakemake.output.formatted_casoff_out)
	offtarget_scores_out = str(snakemake.output.offtarget_scores_out)
	# === Params ===
	rna_bulge = str(snakemake.params.rna_bulge)
	dna_bulge = str(snakemake.params.dna_bulge)
	maximum_mismatches = str(snakemake.params.max_mismatch)
	PU = str(snakemake.params.casoff_accelerator)

	# == Compile Casoff parameters
	casoff_params = (maximum_mismatches, rna_bulge, dna_bulge, PU)

	# == Import per-editor report
	# concatenation using a list of DataFrames
	df_list = []
	for report_path in guides_report_per_editor_list:
		loop_df = pd.read_csv(report_path)
		df_list.append(loop_df)

	# Concatenate all DataFrames at once
	concat_guides_report = pd.concat(df_list)

	# == Import scores dictionary
	concat_offtarget_scores_dict = {}
	for offtarget_scores in offtarget_scores_list:
		with open(offtarget_scores, 'rb') as pickle_handle:
			offtarget_scores_dict = pickle.load(pickle_handle)
			# Update the main dictionary in-place without re-assigning
			concat_offtarget_scores_dict.update(offtarget_scores_dict)

	write_out_res(concat_offtarget_scores_dict,
				  concat_guides_report,
				  casoff_params,
				  offtarget_scores_out,
				  formatted_casoff_out)


if __name__ == "__main__":
	main()

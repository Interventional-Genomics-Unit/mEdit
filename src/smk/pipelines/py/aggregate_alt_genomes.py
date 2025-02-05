# == Native Modules
from decimal import Decimal
import logging
# == Installed Modules
import pandas as pd
import numpy as np
# == Project Modules


def compile_tables(dataframes_list, float_cols):
	# Initialize an empty DataFrame to store aggregated results
	all_data = pd.DataFrame()
	# Loop through dataframes
	for file_path in dataframes_list:
		df = pd.read_csv(file_path)
		# Replace '-' with NaN and let pandas guess data types
		for col_name in float_cols:
			try:
				df[col_name] = df[col_name].replace('-', np.nan).astype(float)
			except KeyError:
				continue
		all_data = pd.concat([all_data, df], ignore_index=True)
	return all_data

# def consolidate_values(series):
# 	"""Consolidate concordant values and store divergent values as a comma-separated string."""
# 	if series.dtype == 'float':
# 		# Convert float values to Decimal for precise comparison
# 		series = series.apply(lambda x: Decimal(str(x)) if not pd.isna(x) else x)
# 	unique_values = series.dropna().unique()
# 	if len(unique_values) == 1:
# 		return unique_values[0]  # Concordant value
# 	return ','.join(map(str, unique_values))  # Divergent values

def consolidate_values(series):
	"""Consolidate concordant values and store divergent values as a comma-separated string."""
	processed = []
	for val in series:
		if pd.isna(val):
			continue  # Skip NaN values
		# Convert numeric values to Decimal for precise comparison
		if isinstance(val, (int, float)):
			# Handle floats and integers as Decimal
			dec_val = Decimal(str(val)).normalize()
			# Remove trailing zeros for consistency (e.g., 5.0 becomes 5)
			dec_val = dec_val.quantize(Decimal(1)) if dec_val == dec_val.to_integral() else dec_val
			processed.append(dec_val)
		else:
			# Treat non-numeric values as strings
			processed.append(str(val))

	unique_values = []
	seen = set()
	for val in processed:
		# Ensure uniqueness based on string representation
		str_val = str(val)
		if str_val not in seen:
			seen.add(str_val)
			unique_values.append(val)

	if len(unique_values) == 0:
		return None  # All values were NaN
	elif len(unique_values) == 1:
		return unique_values[0]  # Return the original type (Decimal, str, etc.)
	else:
		# Join divergent values as strings
		return ','.join(map(str, unique_values))


def aggregate_tables(dataframe):
	"""Aggregate a DataFrame, grouping by Guide_ID and QueryTerm (or QueryTerm alone)."""
	if dataframe.empty:
		return dataframe

	try:
		# Group by Guide_ID and QueryTerm
		aggregated = dataframe.groupby(['Guide_ID', 'QueryTerm'], as_index=False).agg(
			lambda x: consolidate_values(x) if x.name not in ['Guide_ID', 'QueryTerm'] else x.iloc[0]
		)
	except KeyError:
		# Fallback to grouping by QueryTerm alone
		aggregated = dataframe.groupby('QueryTerm', as_index=False).agg(
			lambda x: consolidate_values(x) if x.name != 'QueryTerm' else x.iloc[0]
		)

	return aggregated


# def aggregate_tables(dataframe):
# 	# Handle empty DataFrame
# 	if dataframe.empty:
# 		return dataframe
# 	try:
# 		aggregated_data = dataframe.groupby(['Guide_ID', 'QueryTerm']).agg(
# 			lambda x: consolidate_values(x) if x.name != 'Guide_ID' else x.iloc[0]).reset_index()
# 	except KeyError:
# 		aggregated_data = dataframe.groupby('QueryTerm').agg(
# 			lambda x: consolidate_values(x) if x.name != 'QueryTerm' else x.iloc[0]).reset_index()
# 	return aggregated_data


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	diffguides_out_list = list(snakemake.input.diff_guides)
	altvar_out_list = list(snakemake.input.alt_var)
	# === Outputs ===
	aggr_diff_guides = str(snakemake.output.aggr_diff_guides)
	aggr_alt_var = str(snakemake.output.aggr_alt_var)
	#   == Log File ==
	logfile_path = str(snakemake.output.logfile_path)
	# === Params ===
	float_cols = list(snakemake.params.float_cols)

	# #DEBUG
	# diffguides_out_list = ['/Users/bellieny/projects/mEdit/dump/guides_report_HG02886/0_Guide_differences.csv',
	# 					   '/Users/bellieny/projects/mEdit/dump/guides_report_HG03453/0_Guide_differences.csv',
	# 					   '/Users/bellieny/projects/mEdit/dump/guides_report_HG02622/0_Guide_differences.csv']
	# altvar_out_list = ['/Users/bellieny/projects/mEdit/dump/guides_report_HG02886/0_Guide_differences.csv',
	# 					   '/Users/bellieny/projects/mEdit/dump/guides_report_HG03453/0_Guide_differences.csv',
	# 					   '/Users/bellieny/projects/mEdit/dump/guides_report_HG02622/0_Guide_differences.csv']
	# float_cols = ['Alt CBE Score', 'Alt ABE score',
	# 			  'Alt Azimuth Score', 'Alt DeepCas9 Score',
	# 			  'Alt DeepCpf1 Score', 'Alt OOF Score']

	# === Log Process Initialization
	# Configure the logging system
	logging.basicConfig(
		level=logging.DEBUG,  # Set the minimum log level (DEBUG logs everything)
		format='%(asctime)s [%(levelname)s] %(message)s',  # Define log format
		handlers=[
			logging.FileHandler(logfile_path),  # Log to file
		]
	)

	logging.info('=== INITIALIZING AGGREGATION ROUTINE FOR ALTERNATIVE GENOMES REPORTS ===')

	# == Compile guide differences and alternative variants tables from different alternative genomes
	diffguides_full_df = compile_tables(diffguides_out_list, float_cols)
	altvar_full_df = compile_tables(altvar_out_list, float_cols)

	# === Apply the consolidation logic for all columns except QueryTerm
	aggregated_diff_data = aggregate_tables(diffguides_full_df)
	aggregated_altvar_data = aggregate_tables(altvar_full_df)

	# === Save to CSV
	aggregated_diff_data.to_csv(aggr_diff_guides, index=False)
	aggregated_altvar_data.to_csv(aggr_alt_var, index=False)

	logging.info('=== AGGREGATION ROUTINE FOR ALTERNATIVE GENOMES REPORTS FINALIZED ===')


if __name__ == "__main__":
	main()

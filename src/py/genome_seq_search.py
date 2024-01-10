# Native Modules
import subprocess
# import regex as re
# import sys
import os
from datetime import date  # datetime, date
# Installed Modules
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from Bio.Blast.Applications import NcbimakeblastdbCommandline as makedb
from Bio.Blast import NCBIXML
# from Bio import motifs


def set_export(outdir):
	# Create outdir if inexistent
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	return outdir


def blast_id_parsing(id_string, sep_dict):
	loop_id = id_string
	for level in sep_dict:
		try:
			loop_id = loop_id.split(sep_dict[level][0])[sep_dict[level][1]]
		except IndexError:
			continue
	return loop_id


def hit_report_limit(query_to_hit_dict, nkeep):
	filtered_dict = {}
	hit_id_list = []
	# Within each blast query, keep a maximum of nkeep hits based on lowest evaue
	for query in query_to_hit_dict:
		series = pd.DataFrame.from_dict(query_to_hit_dict[query], orient='index')
		fseries = series['e-value'].nsmallest(nkeep, keep='first')
		for hit in fseries.index.tolist():
			# Keep the same structure of the input dictionary
			filtered_dict.setdefault(query, {}).setdefault(hit, query_to_hit_dict[query][hit])
			hit_id_list.append(hit)
	return filtered_dict, hit_id_list


def blast_result_parser(xml_temp_path, eval_threshold, sep_instructions, nkeep):
	report_dict = {}
	for record in NCBIXML.parse(open(xml_temp_path)):
		if record.alignments:
			for align in record.alignments:
				for hsp in align.hsps:
					if hsp.expect < eval_threshold:
						id = blast_id_parsing(align.hit_def, sep_instructions)
						report_dict.setdefault(record.query, {}).setdefault(id, {
							"blastp_alignment_len": align.length, "e-value": hsp.expect})
	filtered_report_dict, hit_id_list = hit_report_limit(report_dict, nkeep)
	return filtered_report_dict, hit_id_list


def main():
	# SNAKEMAKE IMPORTS
	#   Inputs
	# The aim is to have a multifasta containing all the syntenic
	# regions in the pangenomes to the canonical region in the Hg38 assembly
	# It's crucial that the path where these are stored is isolated for the creation of a BlastDB
	pangenome_multifasta = str(snakemake.input.multifasta)
	# This will be the consensus sequence between patient's VCF and the Hg38 assembly
	consensus_fasta_path = str(snakemake.input.query_fasta)
	#   Outputs
	resultsfolder = set_export(str(snakemake.output))
	#   Params
	db_path = str(snakemake.params.db_path)
	evalue = str(snakemake.params.blast_evalue_thresh)
	threads = int(snakemake.wildcards.threads)
	blastout_path = str(snakemake.params.blastout_path)

	# Internal path variables
	temp = f"{blastout_path}/temp.xml"

	# MakeblastDB
	makedb_cline = makedb(dbtype="prot",
	                      input_file=pangenome_multifasta)
	makedb_cline()
	# Blast run
	print("BlastP Run in progress")
	blast_cline = blastn(query=consensus_fasta_path,
	                     db=db_path,
	                     task='megablast',
	                     num_threads=threads,
	                     evalue=evalue,
	                     out=temp,
	                     outfmt=5)

	print(f"Blast command used:\n {blast_cline}")
	blast_cline()

	# Format parser
	# TODO: Adjust FASTA ID separator and the number of hits to keep
	print("Parsing BlastP output")
	hit_dict, hit_list = blast_result_parser(temp,
	                                         evalue,
	                                         config["id_sep_dict"],
	                                         config["blast_nkeep"])
	# Remove blast XML temp file
	os.remove(temp)


if __name__ == "__main__":
	main()

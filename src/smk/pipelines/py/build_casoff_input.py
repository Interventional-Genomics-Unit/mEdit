# == Native Modules
import pickle
from collections import defaultdict
import shutil
# == Installed Modules
import pandas as pd
# == Project Modules


def cas_offinder_bulge(tmp_processing_filename, casoff_in_path, bulge: bool):
	'''
	 The cas-offinder off-line package contains a bug that doesn't allow bulges
	This script is partially a wrapper for cas-offinder to fix this bug
	 created by...
	https://github.com/hyugel/cas-offinder-bulge
	'''
	# INPUT LEG
	id_dict = {}
	if bulge:
		with open(tmp_processing_filename, 'r') as f:
			chrom_path = f.readline()
			pattern, bulge_dna, bulge_rna = f.readline().strip().split()
			isreversed = False
			for i in range(int(len(pattern) / 2)):
				if pattern[i] == 'N' and pattern[len(pattern) - i - 1] != 'N':
					isreversed = False
					break
				elif pattern[i] != 'N' and pattern[len(pattern) - i - 1] == 'N':
					isreversed = True
					break
			bulge_dna, bulge_rna = int(bulge_dna), int(bulge_rna)
			targets = [line.strip().split() for line in f]
			rnabulge_dic = defaultdict(lambda: [])
			bg_tgts = defaultdict(lambda: set())
			for raw_target, mismatch, gid in targets:
				if isreversed:
					target = raw_target.lstrip('N')
					len_pam = len(raw_target) - len(target)
					bg_tgts['N' * len_pam + target + 'N' * bulge_dna].add(mismatch)
					id_dict['N' * len_pam + target + 'N' * bulge_dna] = gid
					for bulge_size in range(1, bulge_dna+1):
						for i in range(1, len(target)):
							bg_tgt = 'N' * len_pam + target[:i] + 'N' * bulge_size + target[i:] + 'N' * (bulge_dna - bulge_size)
							bg_tgts[bg_tgt].add(mismatch)
							id_dict[bg_tgt] = gid

					for bulge_size in range(1, bulge_rna+1):
						for i in range(1, len(target)-bulge_size):
							bg_tgt = 'N' * len_pam + target[:i] + target[i+bulge_size:] + 'N' * (bulge_dna + bulge_size)
							bg_tgts[bg_tgt].add(mismatch)
							rnabulge_dic[bg_tgt].append((i, int(mismatch), target[i:i+bulge_size]))
							id_dict[bg_tgt] = gid
				else:
					target = raw_target.rstrip('N')
					len_pam = len(raw_target) - len(target)
					bg_tgts['N' * bulge_dna + target + 'N' * len_pam].add(mismatch)
					id_dict['N' * bulge_dna + target + 'N' * len_pam] = gid
					for bulge_size in range(1, bulge_dna+1):
						for i in range(1, len(target)):
							bg_tgt = 'N' * (bulge_dna - bulge_size) + target[:i] + 'N' * bulge_size + target[i:] + 'N' * len_pam
							bg_tgts[bg_tgt].add(mismatch)
							id_dict[bg_tgt] = gid

					for bulge_size in range(1, bulge_rna+1):
						for i in range(1, len(target)-bulge_size):
							bg_tgt = 'N' * (bulge_dna + bulge_size) + target[:i] + target[i+bulge_size:] + 'N' * len_pam
							bg_tgts[bg_tgt].add(mismatch)
							rnabulge_dic[bg_tgt].append( (i, int(mismatch), target[i:i+bulge_size]) )
							id_dict[bg_tgt] = gid
			if isreversed:
				seq_pam = pattern[:len_pam]
			else:
				seq_pam = pattern[-len_pam:]
		with open(casoff_in_path, 'w') as f:
			f.write(chrom_path)
			if isreversed:
				f.write(pattern + bulge_dna*'N' + '\n')
			else:
				f.write(bulge_dna*'N' + pattern + '\n')
			cnt = 0
			for tgt, mismatch in bg_tgts.items():
				f.write(tgt + ' ' + str(max(mismatch)) + ' ' + '\n')
				cnt += 1
		# THIS FILE PATH IS SUPPLIED TO CASOFF-FINDER
	if not bulge:
		nobulge_dict = {}
		with open(tmp_processing_filename, 'r') as inf:
			for line in inf:
				entry = line.strip().split(' ')
				if len(entry) > 2 and len(entry[-1]) > 3:
					seq, mm, gid = entry
					nobulge_dict[seq] = [gid, mm]
					shutil.copy2(tmp_processing_filename, casoff_in_path)
	return seq_pam
	# print("Running Cas-OFFinder (output file: %s)..." % outfn)


def check_bulge(casoff_params):
	if casoff_params[1:3] == (0, 0):
		bulge = False
	else:
		bulge = True
	return bulge


def make_casoffinder_input(tmp_processing_filename, fasta_fname, pam, pamISfirst, guidelen, guides, gnames, casoff_params):
	## create input file for cas-offinder
	mm, RNAbb, DNAbb, PU = casoff_params

	with open(tmp_processing_filename, 'w') as f:
		f.writelines(fasta_fname + "\n")
		line = 'N' * guidelen

		if pamISfirst:
			line = f"{pam}{line} {DNAbb} {RNAbb}" + "\n"
		else:
			line = f"{line}{pam} {DNAbb} {RNAbb}" + "\n"
		f.writelines(line)
		dpam = 'N' * len(pam)
		for grna, gname in zip(guides, gnames):
			if pamISfirst:
				line = f"{dpam}{grna} {mm} {gname}" + "\n"
			else:
				line = f"{grna}{dpam} {mm} {gname}" + "\n"
			f.writelines(line)


def main():
	'''
	  ### For Daniel to snakemake <---------------

	  input paths for this script:
		  -resultsfolder -- results output folder
		  -guide_search_params -- search paramters used in fetchguides
		  -guide_tab_fname -- original guide table output from FetchGuides OR ALT process_genomes files
		  -fasta_fname -- Hg38 fasta or if using alternative consensus genome
		  -(maybe?) casoffinder path

	  input variables for the script:
		  -genome_name -- name of fasta/consensus we are searching
		  -guides_src_name -- name of the guides source genome ex. HG38 or HG02257

	  *possibly allow for changes in the cas-offinder parameters
	  see bottom of page
	  '''
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	guides_report_per_editor_path = str(snakemake.input.guides_per_editor_path)
	guide_search_params = str(snakemake.input.guide_search_params)
	assembly_reference_path = str(snakemake.input.assembly_path)
	# snv_site_info = str(snakemake.input.snv_site_info)
	# annote_path = str(snakemake.params.annote_path)
	# === Outputs ===
	casoff_input_path = str(snakemake.output.casoff_input)
	seq_pam_path = str(snakemake.output.seq_pam_path)
	# === Params ===
	tmp_processing_casoff_path = str(snakemake.params.tmp_processing_casoff)
	rna_bulge = str(snakemake.params.rna_bulge)
	dna_bulge = str(snakemake.params.dna_bulge)
	maximum_mismatches = str(snakemake.params.max_mismatch)
	PU = str(snakemake.params.casoff_accelerator)
	# === Wildcards ===
	editing_tool = str(snakemake.wildcards.editing_tool)
	# fastaref_name = str(snakemake.wildcards.sequence_id)
	# resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"

	# paths = listdir(resultsfolder)

	# === Guide search params ===
	search_params = pickle.load(open(guide_search_params, 'rb'))
	guides_report_per_editor = pickle.load(open(guides_report_per_editor_path, 'rb'))
	# search_params = pickle.load(open(resultsfolder + "guide_search_params.pkl", 'rb'))

	# hg38 or consensus sequence
	fasta_fname = assembly_reference_path
	# genome_name = fastaref_name
	# fasta_fname = '/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa'
	# genome_name = 'hg38'

	# hg38 guides found (but could be {alt_genome}_differences.csv
	# guide_tab_fname = guides_report_per_editor
	# guides_src_name = guideref_name
	# guide_tab_fname = resultsfolder + 'Guides_found.csv'
	# guide_src_name = 'hg38'

	### Daniel---> Pycharm is not find subprocess.Popen(casoffinder...) without an absolute path. so I'm adding this
	# but I don't think its needed in the final version
	# cas_off_expath = '/home/thudson/miniconda3/envs/edit/bin/cas-offinder'

	# defaults - we may allow users to change these cas-offinder settings?
	# according to Gorodkin et al. and Lin et al.  DNA bulges are even more tolerated than mismatches alone
	# https://www.nature.com/articles/s41467-022-30515-0
	# RNAbb = 0  # RNA bulge, a deletion in the off-target
	# DNAbb = 1  # DNA bulge, an insertion in the off-target
	# mm = 3  # max allowable mismatch
	# PU = 'C'  # G = GPU C = CPU A = Accelerators -- I don't really know which should be default?
	# casoff_params = (3, 0, 0, "C")
	pam, pamISfirst, guidelen = search_params[editing_tool][0:3]
	guides, gnames = list(guides_report_per_editor.gRNA), list(guides_report_per_editor.Guide_ID)

	#
	casoff_params = (maximum_mismatches, rna_bulge, dna_bulge, PU)
	bulge_check = check_bulge(casoff_params)

	make_casoffinder_input(tmp_processing_casoff_path,
						   fasta_fname,
						   pam,
						   pamISfirst,
						   guidelen,
						   guides,
						   gnames,
						   casoff_params)

	# cas_offinder_bulge FUNCTION STARTS HERE
	seq_pam = cas_offinder_bulge(tmp_processing_casoff_path, casoff_input_path, bulge_check)
	with open(seq_pam_path, 'wb') as seq_pam_handle:
		pickle.dump(seq_pam, seq_pam_handle)


if __name__ == "__main__":
	main()

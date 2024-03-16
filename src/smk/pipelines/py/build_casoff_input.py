# == Native Modules
import pickle
from collections import defaultdict
import shutil
# == Installed Modules
import pandas as pd
# == Project Modules


def parse_bulge(tmp_processing_filename):
	isreversed = False
	id_dict = {}
	with open(tmp_processing_filename, 'r') as f:
		chrom_path = f.readline()
		pattern, bulge_dna, bulge_rna = f.readline().strip().split()
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
				for bulge_size in range(1, bulge_dna + 1):
					for i in range(1, len(target)):
						bg_tgt = 'N' * len_pam + target[:i] + 'N' * bulge_size + target[i:] + 'N' * (
									bulge_dna - bulge_size)
						bg_tgts[bg_tgt].add(mismatch)
						id_dict[bg_tgt] = gid

				for bulge_size in range(1, bulge_rna + 1):
					for i in range(1, len(target) - bulge_size):
						bg_tgt = 'N' * len_pam + target[:i] + target[i + bulge_size:] + 'N' * (bulge_dna + bulge_size)
						bg_tgts[bg_tgt].add(mismatch)
						rnabulge_dic[bg_tgt].append((i, int(mismatch), target[i:i + bulge_size]))
						id_dict[bg_tgt] = gid
			else:
				target = raw_target.rstrip('N')
				len_pam = len(raw_target) - len(target)
				bg_tgts['N' * bulge_dna + target + 'N' * len_pam].add(mismatch)
				id_dict['N' * bulge_dna + target + 'N' * len_pam] = gid
				for bulge_size in range(1, bulge_dna + 1):
					for i in range(1, len(target)):
						bg_tgt = 'N' * (bulge_dna - bulge_size) + target[:i] + 'N' * bulge_size + target[
																								  i:] + 'N' * len_pam
						bg_tgts[bg_tgt].add(mismatch)
						id_dict[bg_tgt] = gid

				for bulge_size in range(1, bulge_rna + 1):
					for i in range(1, len(target) - bulge_size):
						bg_tgt = 'N' * (bulge_dna + bulge_size) + target[:i] + target[i + bulge_size:] + 'N' * len_pam
						bg_tgts[bg_tgt].add(mismatch)
						rnabulge_dic[bg_tgt].append((i, int(mismatch), target[i:i + bulge_size]))
						id_dict[bg_tgt] = gid
		if isreversed:
			seq_pam = pattern[:len_pam]
		else:
			seq_pam = pattern[-len_pam:]
	return isreversed, chrom_path, seq_pam, rnabulge_dic, id_dict, len_pam, pattern, bulge_dna, bg_tgts


def write_casoff_input(tmp_processing_filename, casoff_in_path, bulge_check: bool):
	'''
	 The cas-offinder off-line package contains a bug that doesn't allow bulges
	This script is partially a wrapper for cas-offinder to fix this bug
	 created by...
	https://github.com/hyugel/cas-offinder-bulge
	'''
	if bulge_check:
		(isreversed,
		 chrom_path,
		 seq_pam,
		 rnabulge_dic,
		 id_dict,
		 len_pam,
		 pattern,
		 bulge_dna,
		 bg_tgts) = parse_bulge(tmp_processing_filename)

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
	if not bulge_check:
		shutil.copy2(tmp_processing_filename, casoff_in_path)


def check_bulge(casoff_params):
	if casoff_params[1:3] == (0, 0):
		bulge = False
	else:
		bulge = True
	return bulge


def parse_casoff_support(tmp_processing_filename, fasta_fname, pam, pamISfirst, guidelen, guides, gnames, casoff_params):
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
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	guides_report_per_editor_path = str(snakemake.input.guides_per_editor_path)
	guide_search_params = str(snakemake.input.guide_search_params)
	assembly_reference_path = str(snakemake.input.decompressed_assembly_symlink)
	# snv_site_info = str(snakemake.input.snv_site_info)
	# annote_path = str(snakemake.params.annote_path)
	# === Outputs ===
	casoff_support_path = str(snakemake.output.casoff_support)
	casoff_input_path = str(snakemake.output.casoff_input)
	# === Params ===
	rna_bulge = str(snakemake.params.rna_bulge)
	dna_bulge = str(snakemake.params.dna_bulge)
	maximum_mismatches = str(snakemake.params.max_mismatch)
	PU = str(snakemake.params.casoff_accelerator)
	# === Wildcards ===
	editing_tool = str(snakemake.wildcards.editing_tool)

	# === Guide search params ===
	search_params = pickle.load(open(guide_search_params, 'rb'))
	guides_report_per_editor = pickle.load(open(guides_report_per_editor_path, 'rb'))

	pam, pamISfirst, guidelen = search_params[editing_tool][0:3]
	guides, gnames = list(guides_report_per_editor.gRNA), list(guides_report_per_editor.Guide_ID)

	# Check bulge based on pre-defined Cas-Offinder params
	casoff_params = (maximum_mismatches, rna_bulge, dna_bulge, PU)
	bulge_check = check_bulge(casoff_params)

	parse_casoff_support(casoff_support_path,
						 assembly_reference_path,
						 pam,
						 pamISfirst,
						 guidelen,
						 guides,
						 gnames,
						 casoff_params)

	# cas_offinder_bulge FUNCTION STARTS HERE
	write_casoff_input(casoff_support_path, casoff_input_path, bulge_check)


if __name__ == "__main__":
	main()

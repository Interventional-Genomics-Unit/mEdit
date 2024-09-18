# == Native Modules
# == Installed Modules
# == Project Modules
import pickle

from build_casoff_input import (check_bulge, parse_bulge)
from scoring import (cfd_score, cfd_spec_score)
from annotate import Transcript


def parse_nobulge(tmp_processing_filename):
	nobulge_dict = {}
	with open(tmp_processing_filename, 'r') as inf:
		for line in inf:
			entry = line.strip().split(' ')
			if len(entry) > 2 and len(entry[-1]) > 3:
				seq, mm, gid = entry
				nobulge_dict[seq] = [gid, mm]
	return nobulge_dict


def write_casoff_output(raw_casoff_out, casoff_support_file, formatted_output, bulge_check):
	with open(raw_casoff_out) as fi, open(formatted_output, 'w') as fo:
		fo.write('Coordinates\tDirection\tGuide_ID\tBulge type\tcrRNA\tDNA\tMismatches\tBulge Size\n') \

		ot_coords = []
		for line in fi:
			entries = line.strip().split('\t')
			ncnt = 0

			if not bulge_check:
				nobulge_dict = parse_nobulge(casoff_support_file)
				gid, mm = nobulge_dict[entries[0]]
				coord = f'{entries[1]}:{entries[2]}-{int(entries[2]) + len(entries[0])}'
				fo.write(f'{coord}\t{entries[4]}\t{gid}\tX\t{entries[0]}\t{entries[3]}\t{entries[5]}\t0\n')
				ot_coords.append(coord)
			if bulge_check:
				(isreversed,
				 chrom_path,
				 seq_pam,
				 rnabulge_dic,
				 id_dict,
				 len_pam,
				 pattern,
				 bulge_dna,
				 bg_tgts) = parse_bulge(casoff_support_file)

				if isreversed:
					for c in entries[0][::-1]:
						if c == 'N':
							ncnt += 1
							break
					if ncnt == 0:
						ncnt = -len(entries[0])
				else:
					for c in entries[0]:
						if c == 'N':
							ncnt += 1
						else:
							break

				if entries[0] in rnabulge_dic:
					gid = id_dict[entries[0]]
					for pos, query_mismatch, seq in rnabulge_dic[entries[0]]:
						if isreversed:
							tgt = (seq_pam + entries[0][len_pam:len_pam + pos] + seq + entries[0][len_pam + pos:-ncnt],
								   entries[3][:len_pam + pos] + '-' * len(seq) + entries[3][len_pam + pos:-ncnt])
						else:
							tgt = (entries[0][ncnt:ncnt + pos] + seq + entries[0][ncnt + pos:-len_pam] + seq_pam,
								   entries[3][ncnt:ncnt + pos] + '-' * len(seq) + entries[3][ncnt + pos:])
						if query_mismatch >= int(entries[5]):
							start = int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (
									isreversed and ncnt > 0 and entries[4] == "-") else 0)
							coord = f'{entries[1]}:{start}-{int(start) + len(tgt[1])}'
							ot_coords.append(coord)
							fo.write(
								f'{coord}\t{entries[4]}\t{gid}\tRNA\t{tgt[0]}\t{tgt[1]}\t{int(entries[5])}\t{len(seq)}\n')

				else:
					gid = id_dict[entries[0]]
					nbulge = 0
					if isreversed:
						for c in entries[0][:-ncnt][len_pam:]:
							if c == 'N':
								nbulge += 1
							elif nbulge != 0:
								break
						tgt = (seq_pam + entries[0][:-ncnt][len_pam:].replace('N', '-'), entries[3][:-ncnt])
					else:
						for c in entries[0][ncnt:][:-len_pam]:
							if c == 'N':
								nbulge += 1
							elif nbulge != 0:
								break
						tgt = (entries[0][ncnt:][:-len_pam].replace('N', '-') + seq_pam, entries[3][ncnt:])
					start = int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (
							isreversed and ncnt > 0 and entries[4] == "-") else 0)
					btype = 'X' if nbulge == 0 else 'DNA'
					coord = f'{entries[1]}:{start}-{int(start) + len(tgt[1])}'
					ot_coords.append(coord)
					fo.write(
						f'{entries[1]}:{start}-{start + len(tgt[1])}\t{entries[4]}\t{gid}\t{btype}\t{tgt[0]}\t{tgt[1]}\t{int(entries[5])}\t{nbulge}\n')

		editor = gid.split('_')[0]
		print(f'{len(ot_coords)} off targets found for {editor}')


def score_ot(crrna,otseq,models_dir):
	score = 0
	if '-' not in crrna[:-3]:
		if '-' not in otseq[:-4]:
			# TODO: NEEDS
			#raise Exception("The function scoring.cfd_score requires a path to the model files")
			score = cfd_score(crrna[:-3].upper(), otseq.upper(),models_dir)
	return score


def annotate_ots(output_filename, annote_path, models_dir):
	'''
	Reads output, Scores each off-target seq and annotates each off_target
	'''

	editor = output_filename.split('_casoff')[0].split('_')[-1]
	coords = []
	scores = []
	spec_scores = {}
	res = open(output_filename, 'r').readlines()
	for line in res[1:]:
		line = line.strip().split('\t')
		coords.append(line[0])
		if line[2] not in spec_scores.keys():
			spec_scores[line[2]] = 0 if editor == 'spCas9' else '.'
		if editor == 'spCas9':
			score = score_ot(line[4], line[5],models_dir)
			if score > 0 and score != 1:
				spec_scores[line[2]] = spec_scores[line[2]] + score
		else:
			score = '.'
		scores.append(score)

	Transcript.load_transcripts(annote_path,coords)

	new_lines = []
	annotate_out = output_filename.replace('_output', '_annotated')
	with open(annotate_out, 'w') as anout:
		anout.write(res[0].strip() + f'\tAnnotation\tScore\n')
		cnt = 0
		for line in res[1:]:
			line = line.strip().split('\t')
			ans = Transcript.transcript(line[0])

			if ans != 'intergenic':
				tid, gname = ans.tx_info()[1:3]
				feature = ans.feature
				x = '|'.join([feature, gname, tid])
				new_line = line + [x, str(scores[cnt])]
				new_lines.append(new_line)

			else:
				x = 'intergenic'
				new_line = line + [x, str(scores[cnt])]
				new_lines.append(new_line)
			cnt += 1
			anout.write('\t'.join(new_line) + '\n')

	if editor == 'spCas9':
		for gid, sum_score in spec_scores.items():
			if sum_score!= 0:
				spec_scores[gid] = cfd_spec_score(sum_score)
		else:
			spec_scores[gid] = "."
	#remove(output_filename)
	return new_lines, spec_scores


def agg_results(lines,mmco):
	'''
	sums the number of off-targets for each guide
	aggregates by mismatch cutoff and bulge size
	'''
	ots_dict = {}
	for line in lines:
		gid, btype, mm, bsize = line[2], line[3], line[6], line[7]
		if gid not in ots_dict.keys():
			ots_dict[gid] = {}
			for i in ['X', 'RNA', 'DNA']:
				for j in range(mmco+1):
					ots_dict[gid][(i,j)] = 0
		ots_dict[gid][(btype, int(mm))] += 1
	return ots_dict


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	raw_casoff_output_path = str(snakemake.input.casoff_out)
	casoff_support_path = str(snakemake.input.casoff_support)
	# === Outputs ===
	offtarget_scores = str(snakemake.output.offtarget_scores)
	formatted_casoff_temp = str(snakemake.output.formatted_casoff_temp)
	# === Params ===
	rna_bulge = str(snakemake.params.rna_bulge)
	dna_bulge = str(snakemake.params.dna_bulge)
	maximum_mismatches = str(snakemake.params.max_mismatch)
	PU = str(snakemake.params.casoff_accelerator)
	annote_path = str(snakemake.params.annote_path)
	models_path = str(snakemake.params.models_path)

	# Check bulge based on pre-defined Cas-Offinder params
	casoff_params = (maximum_mismatches, rna_bulge, dna_bulge, PU)
	bulge_check = check_bulge(casoff_params)

	write_casoff_output(raw_casoff_output_path,
						casoff_support_path,
						formatted_casoff_temp,
						bulge_check)

	# == Process off-target scoring
	new_lines, spec_scores = annotate_ots(formatted_casoff_temp, annote_path, models_path)
	ot_dict = agg_results(new_lines, casoff_params[0])  # sum off-targets

	# == Compile and export scores
	ots = {}
	for k, v in ot_dict.items():
		v['spec_score'] = spec_scores[k]
		ots[k] = v
	with open(offtarget_scores, 'wb') as pickle_handle:
		pickle.dump(ots, pickle_handle)


if __name__ == "__main__":
	main()

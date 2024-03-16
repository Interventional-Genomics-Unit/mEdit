# == Native Modules
# == Installed Modules
# == Project Modules
from build_casoff_input import (check_bulge, parse_bulge)


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


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	raw_casoff_output_path = str(snakemake.input.casoff_out)
	casoff_support_path = str(snakemake.input.casoff_support)
	# === Outputs ===
	formatted_casoff_out = str(snakemake.output.formatted_casoff_out)
	# === Params ===
	rna_bulge = str(snakemake.params.rna_bulge)
	dna_bulge = str(snakemake.params.dna_bulge)
	maximum_mismatches = str(snakemake.params.max_mismatch)
	PU = str(snakemake.params.casoff_accelerator)

	# Check bulge based on pre-defined Cas-Offinder params
	casoff_params = (maximum_mismatches, rna_bulge, dna_bulge, PU)
	bulge_check = check_bulge(casoff_params)

	write_casoff_output(raw_casoff_output_path, casoff_support_path, formatted_casoff_out, bulge_check)


if __name__ == "__main__":
	main()

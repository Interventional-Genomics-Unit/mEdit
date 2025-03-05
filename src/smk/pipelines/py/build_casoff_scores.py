# == Native Modules
# == Installed Modules
import pandas as pd
# == Project Modules
from scoring import cfd_score,load_model_params
from annotate import Transcript

#######
# Reformats guidescan output to include 1) corrected cfd scores 2) distances 3) genome found 4) genomic annotations
# For alternative genomes, eliminates any sites that do not differ from reference genome
# By Taylor H.
######

def fix_cfd(lines,models_dir,score_guides=True):
	'''
	fixes guidescan2 cfd scoring and adds edit distance to each site
	:param lines: lines from guidescan bed file
	:param models_dir: cfd score model weights
	:return:
	'''
	new_lines = []
	new_lines.append(lines[0]+",Distance")
	if score_guides:
		mm_scores, pam_scores = load_model_params('cfd', models_dir)
	for line in lines[1:]:
		line = line.split(",")
		seq1 = line[1][:-3]
		seq2 = line[6][:-3]
		pam = line[6][-3:]
		if score_guides and (int(line[7]) + int(line[8]))== 0:
			score = cfd_score(seq1, seq2, pam,mm_scores, pam_scores)
		else:
			score = -1
		dist = int(line[5]) + int(line[7]) + int(line[8])
		line[-4] = str(score)
		line.append(str(dist))
		new_lines.append(','.join(line))
	return new_lines


def reformat_ref_and_alt(lines,offtarget_genome,genome_type):
	lines_reformatted = []
	header = 'id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity,alt_site_impact,alt_var,alt_genome'
	ref = True if genome_type != 'extended' else False
	lines_reformatted.append(header)
	for line in lines:
		line_split = line.strip().replace("\t",",").replace("\n","").split(",")[3:]
		line_split[6] = line_split[6].replace(".","-")
		line_split[1] = line_split[1].replace(".", "-")
		if ref:
			newline = ",".join(line_split[:10]) + ",na,na,none"
			lines_reformatted.append(newline)
		else:

			variant = line_split[11:]
			#drop lines where variants are not in sequence
			dist_from_variant = (int(line_split[12])) - int(line_split[3])

			if dist_from_variant < len(line_split[6].replace("-",""))-1	 and dist_from_variant >= -1:
				hgvs = f"{variant[0]}:{variant[1]}{variant[3]}>"
				for alt_allele in variant[4:]: #incase multi-allelic
					hgvs += f"{alt_allele}|"
				newline = line_split[:11] + [hgvs[:-1]] +[str(offtarget_genome)]
				lines_reformatted.append(",".join(newline))
	return lines_reformatted


def compile_mutliple_alignments(dist_lines,max_bulge):
	# creates a dict for every sites that has mulitple alignments
	ot_dict = {}
	for line in dist_lines[1:]:
		line = line.strip().split(",")
		pos = int(line[3])
		prefix = line[0]+"_"+line[11]+"_"+line[2] + "_"
		not_found = True
		for i in range(0, max_bulge):

			if f"{prefix}{str(pos+i)}" in ot_dict.keys():
				ot_dict[f"{prefix}{str(pos+i)}"].append(line)
				not_found = False
				break
			elif f"{prefix}{str(pos-i)}" in ot_dict.keys():
				ot_dict[f"{prefix}{str(pos-i)}"].append(line)
				not_found = False
				break
			else:
				pass
		if not_found:
			ot_dict[f"{prefix}{str(pos)}"] = [line]

	return ot_dict


def config_alt_variants(df,find_alt_unique_sites):
	'''
	concatenates variants that span the same site
	'''
	if find_alt_unique_sites:
		cols = [x for x in df.columns if x != "alt_var"]
		new_df = df.groupby(cols, as_index=False).agg({"alt_var": lambda x: "|".join(sorted(set(x)))})

	else:
		cols = [x for x in df.columns if x != "alt_var"] + ['alt_var']
		new_df = df.loc[:,cols]
	return new_df


def de_dup(dist_lines,max_bulge):
	new_lines = []
	header = dist_lines[0].strip().split(",")
	header.append("Alt Alignment [0=No/1=Yes]")
	new_lines.append(header)
	ot_dict = compile_mutliple_alignments(dist_lines,max_bulge)
	for coord,lines in ot_dict.items():
		if len(lines) ==1:
			bestline = lines[0]
			alt_alignment = '0'
		else:
			alt_alignment = '1'
			dist, mm,placement_score = 100,0,0
			bestline = ''

			for line in lines:
				if int(line[-1]) < dist: # if edit distance is lower than other alignments keep
					bestline = line
					dist, mm, placement_score = int(line[-1]), int(line[5]), sum([i for i in range(len(line[6])) if line[6][i].islower() or line[6][i] == "-"])
				elif int(line[-1]) == dist: # if edit distance is equal then keep if mismatch is higher(qalt alignment has more bulges)
					if int(line[5]) > mm:
						bestline = line
						dist, mm, placement_score = int(line[-1]), int(line[5]), sum([i for i in range(len(line[6])) if line[6][i].islower() or line[6][i] == "-"])
					elif int(line[-1]) == mm:  # if edit distance and mm equal then keep if mismatches are further from 3'pam

						if sum([i for i in range(len(line[6])) if line[6][i].islower() or line[6][i] == "-"]) < placement_score:
							bestline = line
							dist, mm, placement_score = int(line[-1]), int(line[5]), sum([i for i in range(len(line[6])) if line[6][i].islower() or line[6][i] == "-"])
					else:
						pass
				else:

					pass
		new_lines.append(bestline+[alt_alignment])
	return new_lines


def add_annotations(df,annote_path):

	coords = df['match_chrm'] + ":" + df['match_position'].astype('str') + "-" + df['match_position'].astype('str')
	Transcript.load_transcripts(annote_path, coords)
	Gene, Feature = [], []
	for snv in coords:
		pos_in_transcript = Transcript.transcript(snv)
		if pos_in_transcript != 'intergenic':
			Gene.append(pos_in_transcript.tx_info()[2])
			Feature.append(pos_in_transcript.feature)
		else:
			Gene.append(".")
			Feature.append('intergenic')

	df['Gene'] = Gene
	df['Feature'] = Feature
	df['match_chrm'] = df['match_chrm'] + ":" + df['match_position'].astype('str') + df['match_strand'].astype('str')

	df = df.iloc[:, [0, 1, 2, 5, 6, 7, 8, 9, 12, 13, 10, 11, 14,15, 16]]
	header = ['Guide_ID', 'On_Target_Sequence', "Match_Coords",'Mismatch',
			  'Match_Sequence', 'RNA_Bulges', 'DNA_Bulges',
			  'CFD_Score', "Distance", 'Multi Alignment [0=No/1=Yes]', "Alt Site Impact",
			  "Alt Genome","Alt Variants",
			 "Feature","Gene"]

	df.columns = header
	df = df.sort_values("Mismatch",ascending = False).sort_values("Distance",ascending = True)
	return df


def reformat_guidescan(guidescan_filtered_bed,
					   formatted_casoff_out,
					   genome_type,
					   offtarget_genome,
					   max_bulge,
					   annote_path,
					   models_dir,
					   editing_tool):

	final_df = pd.DataFrame()
	find_alt_unique_sites = True if genome_type == 'extended' else False
	lines = open(guidescan_filtered_bed,"r").readlines()

	if len(lines) ==1:
		# editing_tool = guidescan_filtered_bed.split('/')[-1].split("_").replace("_guidescan_filtered.bed","").split("_")[-1]
		print(f"No offtargets found for {editing_tool} in {genome_type}")
		print("skipping......")
	else:
		lines_reformatted = reformat_ref_and_alt(lines,offtarget_genome,genome_type)
		scored_lines = fix_cfd(lines_reformatted, models_dir)
		condensed_lines = de_dup(scored_lines,max_bulge)
		df = pd.DataFrame(condensed_lines[1:], columns=condensed_lines[0])
		adjusted_for_variants_df = config_alt_variants(df,find_alt_unique_sites)
		final_df = add_annotations(adjusted_for_variants_df,annote_path)
	final_df.to_csv(formatted_casoff_out,index = False)


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	guidescan_filtered_bed = str(snakemake.input.guidescan_filtered_bed)
	# === Outputs ===
	formatted_casoff_out = str(snakemake.output.formatted_casoff)
	# === Params ===
	models_path = str(snakemake.params.models_path)
	annote_path = str(snakemake.params.annote_path)
	extended_genomes = list(snakemake.params.extended_genomes)
	rna_bulge = int(snakemake.params.rna_bulge)
	dna_bulge = int(snakemake.params.dna_bulge)
	# === Wildcards ===
	offtarget_genome = str(snakemake.wildcards.offtarget_genomes)
	editing_tool = str(snakemake.wildcards.editing_tool)

	genome_type = 'main_ref'
	if offtarget_genome in set(extended_genomes):
		genome_type = 'extended'

	max_bulge = max(rna_bulge, dna_bulge)

	reformat_guidescan(guidescan_filtered_bed,
					   formatted_casoff_out,
					   genome_type,
					   offtarget_genome,
					   max_bulge,
					   annote_path,
					   models_path,
					   editing_tool)


if __name__ == "__main__":
	main()

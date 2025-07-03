# == Native Modules
# == Installed Modules
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
import pickle
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

def revcom(s):
	basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
	letters = list(s[::-1])
	letters = [basecomp[base] for base in letters]
	return ''.join(letters)

def find_codon(pos,tx_seq,strand,rf):
	if strand == '+':
		codon = tx_seq[int(pos - rf): int((pos - rf) + 3)]
	else:
		adj_rf = 2 - rf
		codon = tx_seq[int(pos - adj_rf): int(pos - adj_rf) + 3]
		codon = revcom(codon)
	return codon



def get_AAconversion_type(codon1, codon2 ,aa1, aa2):
	"""
    codon1: codon of Alt allele to be changed by BE
    codon2: codon after changed by BE
    """
	aa_groups = [["G", "A", "V", "L", "I", "M", "F", "Y", "W"],
				 ["S", "Q", "T", "N"],
				 ["C", "G", "P"],
				 ["D", "E"],
				 ["K", "H", "R", "Q"]]
	if aa1 == aa2:
		mtype = 'Synonymous'
		if codon1 == codon2:
			mtype = 'Silent'
	else:
		if codon2 in ['TAA', 'TAG', 'TGA']:
			mtype = 'Nonsense'
		elif ([aa_groups.index(x) for x in aa_groups if str(aa1) in x] ==
			  [aa_groups.index(x) for x in aa_groups if str(aa2) in x]):
			mtype = 'Conservative'
		else:
			mtype = 'Non-conservative'
	return mtype

def get_codon_change(ref_codon,conversion,rf):

	# ref codon
	aa_ref = Seq(ref_codon).translate()
	# new codon
	new_codon = Seq("".join(ref_codon[x] if x != abs(rf) else conversion[1] for x in [0, 1, 2]))
	aa_new = new_codon.translate()

	ctype = get_AAconversion_type(ref_codon, new_codon ,aa_ref, aa_new)

	aa_new = seq3(aa_new, custom_map={"*": "***"})
	aa_ref = seq3(aa_ref, custom_map={"*": "***"})
	aa_change = f"{ctype}({aa_ref}>{aa_new})"
	return aa_change

def add_be_outcomes(tx,fasta,coord,feature, match_sequence, conversion, be_win):

	target_bases = match_sequence.upper().replace("-", "")[be_win[0] - 1:be_win[1] + 1]
	n_bases = target_bases.count(conversion[0])
	convert = 'na'

	if "exon" in feature or "codon" in feature:
		reading_frames = '012' * 20
		if n_bases > 0:
			convert = []
			eid, tid, gname, strand, txstart = tx.tx_info()
			tx_seq = tx.get_tx_seq(fasta).upper()
			for i,base in enumerate(target_bases):
				if base == conversion[0]:
					j = i + (be_win[0] - 1)
					base_pos = int(coord.split("-")[-1]) + j
					t_pos = int(base_pos) - (int(txstart) +1)
					rf = int(reading_frames[reading_frames.index(str(tx.rf)):][j])
					ref_codon = find_codon(t_pos, tx_seq, strand, rf)
					convert.append(get_codon_change(ref_codon,conversion,rf))
		convert = ",".join(convert)

	return n_bases, convert

def add_be_annotations(df,annote_path,guide_params,fasta_path):
	'''
	counting the number of editable bases (be target base) in the be window and
	 determining base change aa consquence

	'''
	pam, pam_is_first,guidelen,be_win = guide_params[0:4]
	conversion = guide_params[-1]

	gene,features,n_editable, codon_change  = [],[],[],[]
	sites_by_chroms ={}
	chroms, starts = [], []

	coords = df['match_chrm'] + ":" + df['match_position'].astype('str') + "-" + df['match_position'].astype('str')
	for k,coord,seq,strand in zip(df['match_chrm'],coords,df.match_sequence,df.match_strand):
		if k not in sites_by_chroms.keys():
			sites_by_chroms[k] ={'coords':[],"match_sequences":[],"strands":[]}
		sites_by_chroms[k]['coords'].append(coord)
		sites_by_chroms[k]['match_sequences'].append(seq)
		sites_by_chroms[k]['strands'].append(strand)
		chroms.append(k)
		starts.append(coord.split(":")[1].split("-")[0])

	for chrom, vals in sites_by_chroms.items():
		fasta_not_loaded = True

		Transcript.load_transcripts(annote_path, vals['coords'])

		for coord,match_sequence,strand in zip(vals['coords'],vals['match_sequences'],vals['strands']):
			tx = Transcript.transcript(coord)
			if tx == 'intergenic':
				gene.append('.')
				feature = tx
			else:
				gene.append(tx.tx_info()[2])
				feature = tx.feature
				if fasta_not_loaded:
					chr_fasta_path = f"{fasta_path}/chr{str(chrom).replace('chr', '')}.pkl"
					with open(chr_fasta_path, 'rb') as pfile:
						fasta = pickle.load(pfile)
					fasta_not_loaded = False

			nbases,convert = add_be_outcomes(tx, fasta,coord, feature, match_sequence, conversion, be_win)
			features.append(feature)
			n_editable.append(nbases)
			codon_change.append(convert)


	annotations_df = pd.DataFrame({'match_chrm':chroms,'match_position':starts,'Gene':gene,'Feature':features, 'N_Bases_Editable':n_editable,'Base_Change_Consequence':codon_change})
	df = df.merge(annotations_df, on=['match_chrm', 'match_position'])

	df['match_chrm'] = df['match_chrm'] + ":" + df['match_position'].astype('str') + df['match_strand'].astype('str')

	df = df.iloc[:, [0, 1, 2, 5, 6, 7, 8, 9,12, 13, 10, 11, 14,15, 16,17,18]]

	header = ['Guide_ID', 'On_Target_Sequence', "Match_Coords",'Mismatch',
			  'Match_Sequence', 'RNA_Bulges', 'DNA_Bulges',
			  'CFD_Score', "Distance", "Alt Variants",'Multi Alignment [0=No/1=Yes]', "Alt Site Impact",
			  "Alt Genome",
			 "Feature","Gene",'N_Base_Editable','Base_Change_Consequence']

	df.columns = header
	df = df.sort_values(["Mismatch","Distance"],ascending = True)
	return df

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
	df.loc[:, ['N_Base_Editable', 'Base_Change_Consequence']] = 'na'
	df['match_chrm'] = df['match_chrm'] + ":" + df['match_position'].astype('str') + df['match_strand'].astype('str')

	df = df.iloc[:, [0, 1, 2, 5, 6, 7, 8, 9,12, 13, 10, 11, 14, 15, 16, 17, 18]]

	header = ['Guide_ID', 'On_Target_Sequence', "Match_Coords", 'Mismatch',
			  'Match_Sequence', 'RNA_Bulges', 'DNA_Bulges',
			  'CFD_Score', "Distance","Alt Variants", 'Multi Alignment [0=No/1=Yes]', "Alt Site Impact",
			  "Alt Genome",
			  "Feature", "Gene", 'N_Base_Editable', 'Base_Change_Consequence']

	df.columns = header
	df = df.sort_values(["Mismatch", "Distance"], ascending=True)
	return df


def reformat_guidescan(guidescan_filtered_bed,
					   formatted_casoff_out,
					   genome_type,
					   offtarget_genome,
					   max_bulge,
					   annote_path,
					   models_dir,
					   editing_tool,
					   guide_params,
					   fasta_path):

	final_df = pd.DataFrame()
	find_alt_unique_sites = True if genome_type == 'extended' else False
	be_flag = True if len(str(guide_params[3])) > 2 else False
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
		if be_flag:
			final_df = add_be_annotations(adjusted_for_variants_df,annote_path,guide_params,fasta_path)
		else:
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
	guide_params = list(snakemake.params.guide_params)
	# === Wildcards ===
	offtarget_genome = str(snakemake.wildcards.offtarget_genomes)
	editing_tool = str(snakemake.wildcards.editing_tool)
	fasta_path = str(snakemake.input.assembly_dir_path) #pickled fasta reference file (even if its extended/alt genome)


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
					   editing_tool,
					   guide_params,
					   fasta_path)


if __name__ == "__main__":
	main()

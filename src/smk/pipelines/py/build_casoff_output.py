# == Native Modules
import numpy as np
import pickle
import re
# == Installed Modules
import pandas as pd
import copy
# == Project Modules

####
# creates aggregated genomic tables
# A) guide cumilative CFD, mismatches and total sites
# cols = Query,Guide_ID,Coords


class Offtarget_Aggregator():
	def __init__(self,
				 guide_id,
				 offtarget_genome,
				 genome_type,
				 thresholds,
				 on_target_coords):
		
		'''
		keeps aggregate scores and total site counts for a specific guide_id and genome searched
		:param guide_id: unique guide_id
		:param offtarget_genome: 
		:param genome_type: extended or main-ref
		:param thresholds: max_mismatches,max_dna_bulge,max_rna_bulge,max both bulges
		:param on_target_coords: ex: 16:1060000
		'''
		# inputs
		self.gid = guide_id
		self.thresholds = thresholds
		self.offtarget_genome = offtarget_genome
		self.genome_type = genome_type
		self.on_target_pos = int(on_target_coords.split(':')[1].split("-")[0])
		self.chrom = on_target_coords.split(':')[0]

		#outputs
		empty_dict = self.init_site_count_dict()
		self.site_counts = empty_dict
		self.cfd_sums = {}
		for k in empty_dict.keys():
			self.cfd_sums[k] = -1

		self.feature_counts = self.init_feature_count_dict()
		self.ref_counts_set = False

	def __str__(self):
		return f"OfftargetAggregator|{self.gid}_{self.offtarget_genome}_{sum(list(self.site_counts.values()))}"

	def __repr__(self):
		return f"OfftargetAggregator|{self.gid}_{self.offtarget_genome}_{sum(list(self.site_counts.values()))}"

	def set_ref_counts(self,ots_counts,feature_counts,cfd_sums):
		if self.ref_counts_set:
			print('ref counts have already been set')
			pass
		else:
			self.cfd_sums= copy.deepcopy(cfd_sums)
			self.site_counts = copy.deepcopy(ots_counts)
			self.feature_counts = copy.deepcopy(feature_counts)
			self.ref_counts_set = True

	def format_gene_feature(self,feature):
		#'stop/start codon, exon, utr, flanking,intron,ncRNA,intergenic'
		feature = feature.replace('non-coding ',"nc").split("-")[0].split(" ")[0]
		if "codon" in feature:
			feature = 'stop/start_codon'
		elif feature[1:]== 'utr':
			feature = feature[1:]
		elif feature == 'non':
			feature = 'intergenic'
		else:
			pass
		return feature


	def add_site(self,site_thresholds,site_coords,cfd,impact,feature):
		site_pos = int(site_coords.split(":")[1][:-1])
		site_chrom = site_coords.replace('chr', '').split(":")[0]
		feature_point,point,cfd_add  = 0,0,0
		feature = self.format_gene_feature(feature)
		
		if site_chrom == self.chrom and self.on_target_pos >= site_pos - 10 and self.on_target_pos <= site_pos + 10:
				feature_point, point, cfd_add = 0, 0, 0

		elif impact == 'removed': # for alt genomes where a site was no longer present, reduce tallies
			point =  -1
			feature_point = -1
			if cfd > 0:
				cfd_add = -1*cfd
		else: # for ref genome and alt genomes where a sitewas added, tally
			point = 1
			feature_point = +1
			if cfd > 0:
				cfd_add =cfd

		self.site_counts[site_thresholds] += point
		self.feature_counts[site_thresholds][feature] += feature_point
		self.cfd_sums[site_thresholds]+=cfd_add

	def init_feature_count_dict(self):
		# creates empty dict for tallying gene annotation features for each site at mm and bulge thresholds
		site_counts = {}
		max_mm, max_db, max_rb, max_both = self.thresholds
		for m in range(0, max_mm + 1):
			for i in range(0, max_db + 1):
				for j in range(0, max_rb + 1):
					site_counts[str(m) + str(i) + str(j)] = {'stop/start_codon':0,
															 'exon':0,
															 'utr':0,
															 "flanking":0,
															 "intron":0,
															 "ncRNA":0,
															 "intergenic":0}
		return site_counts

	def init_site_count_dict(self):
		# creates empty dict for tallying each site at mm and bulge thresholds
		site_counts = {}
		max_mm, max_db, max_rb, max_both = self.thresholds
		for m in range(0, max_mm + 1):
			for i in range(0, max_db + 1):
				for j in range(0, max_rb + 1):
					site_counts[str(m) + str(i) + str(j)] = 0
		return site_counts

	def get_cfd_aggregate_scores(self):
		# reports/formats of CFDscore post summation
		cfd_aggregated = {}
		for site_threshold, cfd_sum in self.cfd_sums.items():
			cfd_aggregated[site_threshold] = str(round(1 / (1 + cfd_sum ), 4)) if cfd_sum != -1 else -1
		return cfd_aggregated

	def as_rows(self):
		#output fed into offtarget library
		# header = ['Guide_ID','Offtarget Genome', 'Number of Mismatches','Number of RNA Bulges','Number of DNA Bulges',
		# 'Total Site Count','CFD/Guidescan Spec Score','Stop/Start Codon Site Count', 'Exon Site Count', 'Utr Site Count',
		# 'Flanking Site Count','Intron Site Count','ncRNA Site Count','Intergenic Site Count']
		stats= []
		cfd_aggregated = self.get_cfd_aggregate_scores()
		for site_threshold in self.site_counts.keys():
			count = self.site_counts[site_threshold]
			score = cfd_aggregated[site_threshold]
			f_counts = list(self.feature_counts[site_threshold].values())

			stats.append([self.gid,self.offtarget_genome,(site_threshold[0]),int(site_threshold[1]),int(site_threshold[2]),count,score]+f_counts)
		return stats


class OffTargets_Library:
	def __init__(self,
				 thresholds,coords_per_guide,offtarget_genomes):
		'''
		Extractes information from off_target_found.csv and compiles into a library/dictionary of Offtarget_Agresggator instances
		'''
		#inputs
		self.thresholds = thresholds
		self.offtarget_genomes = offtarget_genomes
		self.coords_per_guide = coords_per_guide
		self.gids_searched = list(coords_per_guide.keys())
		self.ref_genome = [k for k in offtarget_genomes.keys() if offtarget_genomes[k] != 'extended'][0]
		#outputs
		self.gid_stats = self.populate_empty_library()
		# {'spCas9_0': {'hg38_GCA_000001405.15': Offtarget_Compilier(), 'HG02557' :Offtarget_Compilier()....}

	def add_data_from_file(self,offtarget_csv,offtarget_genome,genome_type):
		# 'Guide_ID,Match_Coords,Mismatch,RNA_Bulges,DNA_Bulges,Alt Site Impact,Feature'

		data = np.genfromtxt(offtarget_csv, delimiter=',',usecols=(0,2,3,5,6,10,14),dtype=str,  skip_header=1)
		cfd_scores = np.genfromtxt(offtarget_csv, delimiter=',', usecols=(7), skip_header=1)

		# Safety measure: Handle empty file case
		if data.size == 0:
			return

		# Ensure data is at least 2D
		if data.ndim == 1:
			data = np.expand_dims(data, axis=0)

		site_thresholds = np.char.add(np.char.add(data[:, 2], data[:, 3]), data[:, 4])
		nrows = data.shape[0]

		if nrows != 0:
			# For alternative genomes set intial counts to the ref sites which do not differ.
			# Subtract and add based on the sites that do differ
			if genome_type == "extended":

				for guide_id in set(data[:,0]):
					# SAFETY MEASURE -> Some Guide_ids from the "Offtargets_found.csv" were absent in guidescan Inputs

					aggregator =self.gid_stats[guide_id][offtarget_genome]

					if aggregator.ref_counts_set:
						pass
					else:
						ref_aggregator = self.gid_stats[guide_id][self.ref_genome]
						ref_counts = copy.deepcopy(ref_aggregator.site_counts)
						ref_cfd_sums = copy.deepcopy(ref_aggregator.cfd_sums)
						ref_features = copy.deepcopy(ref_aggregator.feature_counts)
						aggregator.set_ref_counts(ref_counts, ref_features,ref_cfd_sums)

			for i in range(nrows):
				# SAFETY MEASURE -> Some Guide_ids from the "Offtargets_found.csv" were absent in guidescan Inputs
				try:
					aggregator = self.gid_stats[data[i,0]][offtarget_genome]
					aggregator.add_site(site_thresholds[i],data[i,1],cfd_scores[i],data[i,5],data[i,6])
					self.gid_stats[data[i, 0]][offtarget_genome] = aggregator

				except KeyError:
					continue

	def populate_empty_library(self):
		gid_stats = {}

		for guide_id in self.gids_searched:
			gid_stats[guide_id] = {}
			on_target_coords = self.coords_per_guide[guide_id]
			for offtarget_genome, genome_type in self.offtarget_genomes.items():
				gid_stats[guide_id][offtarget_genome] =Offtarget_Aggregator(guide_id,
																		   offtarget_genome,
																		   genome_type,
																		   self.thresholds,
																			on_target_coords)
		return gid_stats

	def generate_rows(self):
		header = ['Guide_ID','Offtarget Genome', 'Number of Mismatches','Number of RNA Bulges','Number of DNA Bulges',
		 'Total Site Count','CFD/Guidescan Spec Score','Stop/Start Codon Site Count', 'Exon Site Count', 'UTR Site Count',
		 'Flanking Site Count','Intron Site Count','ncRNA Site Count','Intergenic Site Count']
		rows = [header]
		for gid in self.gid_stats.keys():
			for aggregator in self.gid_stats[gid].values():
				rows +=aggregator.as_rows()
		return rows


def extract_genome_name(filepath: str, pattern: str):
	return re.search(pattern, filepath).group(0)


def compile_genome_types(list_of_files, extented_genomes, reference_genome):
	genome_type_dict = {}
	genome_type_dict.setdefault(reference_genome, 'main_ref')
	for extended_genome in extented_genomes:
		for filepath in list_of_files:
			try:
				genome_name = extract_genome_name(filepath, extended_genome)
				genome_type_dict.setdefault(genome_name, 'extended')
			except AttributeError:
				continue
	return genome_type_dict


def create_count_table(ots_dict):
	df = pd.DataFrame(ots_dict)
	spec_score_table = df.T.sort_index()['cfd']
	df = df.drop(index='cfd').fillna(0).astype('int')
	df.index = pd.MultiIndex.from_tuples(list(df.index),names=['BulgeType', 'Number of Mismatches', 'Number of Bulges'])
	expanded_table = df.reset_index()
	return expanded_table, spec_score_table


def create_summary_report(expanded_offtarget_report,guides_report):
	df = expanded_offtarget_report.copy()

	# find aggrgate scores of all sites fourn in a genome
	df['CFD/Guidescan Spec Score'] = df['CFD/Guidescan Spec Score'].astype('float').apply(lambda x: (1/x - 1) if x > 0 else x)
	cfd_df = df[['Guide_ID', 'Offtarget Genome', 'CFD/Guidescan Spec Score']].groupby(
		['Guide_ID', 'Offtarget Genome']).apply(lambda x: 1 / (1 + x[x > 0].sum()))
	scores = cfd_df['CFD/Guidescan Spec Score'].to_list()

	# find summation of annotation type in a genome (only for mismatches, not bulges)
	df = df.loc[(df['Number of RNA Bulges']==0) & (df['Number of DNA Bulges']==0),:]
	df = df.drop(columns=['Number of RNA Bulges', 'Number of DNA Bulges','CFD/Guidescan Spec Score'])
	df = df.groupby(['Guide_ID', 'Offtarget Genome', 'Number of Mismatches']).sum()
	annotation_totals = df.iloc[:, 1:].groupby(level=[0, 1]).sum()

	# find cumulative sites for each mismatch bracket (only for mismatches, not bulges)
	df = df['Total Site Count'].astype('str').groupby(level=[0,1]).apply('|'.join).to_frame()

	df['CFD/Guidescan Spec Score'] = scores
	df = df.join(annotation_totals).reset_index()
	summary_report = guides_report[['QueryTerm', 'GeneName', 'Editor', 'Guide_ID', 'Alt Guide Impact',
       'Alt Genome', 'Coordinates', 'Strand', 'gRNA', 'Pam']].join(df.set_index('Guide_ID'),on='Guide_ID')

	return summary_report


def compile_per_editor_input(guides_per_editor_list):

	guides_report = pd.DataFrame()

	for guides_per_editor in guides_per_editor_list:
		guides_report_per_editor = pickle.load(open(guides_per_editor, 'rb'))
		guides_report = pd.concat([guides_report, guides_report_per_editor])

	coords_per_guide = guides_report[['Guide_ID', 'Coordinates']].set_index('Guide_ID').to_dict()['Coordinates']
	return coords_per_guide, guides_report


def aggregate_guidescan_results(guides_per_editor_list, formatted_casoff_list,
								offtarget_genomes, thresholds,
								offtarget_summary_file, offtarget_expanded_summary_file):

	coords_per_guide, guides_report = compile_per_editor_input(guides_per_editor_list)

	offtarget_lib = OffTargets_Library(thresholds, coords_per_guide, offtarget_genomes)

	for formatted_casoff in formatted_casoff_list:
		offtarget_genome = formatted_casoff.split("/")[-2]
		genome_type = offtarget_genomes[offtarget_genome]
		offtarget_lib.add_data_from_file(formatted_casoff, offtarget_genome, genome_type)

	rows = offtarget_lib.generate_rows()
	expanded_offtarget_report = pd.DataFrame(rows[1:],columns=rows[0])

	expanded_offtarget_report.to_csv(offtarget_expanded_summary_file,index=False)

	summary_report = create_summary_report(expanded_offtarget_report,guides_report)
	summary_report.to_csv(offtarget_summary_file,index = False)


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	guides_per_editor_list = list(snakemake.input.guides_per_editor_list)
	formatted_casoff_list = list(snakemake.input.formatted_casoff_list)
	# === Outputs ===
	offtarget_summary_file = str(snakemake.output.offtarget_summary)
	offtarget_summary_expanded_file = str(snakemake.output.off_target_summary_expanded)
	# === Params ===
	extended_genomes = list(snakemake.params.extended_genomes)
	maximum_mismatches = int(snakemake.params.max_mismatch)
	rna_bulge = int(snakemake.params.rna_bulge)
	dna_bulge = int(snakemake.params.dna_bulge)
	# === Wildcards ===
	reference_genome = str(snakemake.wildcards.reference_id)

	# === Compile genome types of all the genomes ('extended' and 'main_ref')
	offtarget_genomes = compile_genome_types(formatted_casoff_list, extended_genomes, reference_genome)

	thresholds = (maximum_mismatches, dna_bulge, rna_bulge, dna_bulge+rna_bulge)

	aggregate_guidescan_results(guides_per_editor_list,
								formatted_casoff_list,
								offtarget_genomes,
								thresholds,
								offtarget_summary_file,
								offtarget_summary_expanded_file
								)


if __name__ == "__main__":
	main()




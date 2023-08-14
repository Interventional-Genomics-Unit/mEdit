import os
import pandas as pd
import bbi
import regex as re
import argparse
import sys
import subprocess


# hgvs_id = 'NM_000152.5:c.271G>A' ,win_size=[4,8])
query = str(sys.argv[1])
qtype = str(sys.argv[2])

class DataHandler:

	'''
	To call tables and handle clinvar updates
	'''

	def __init__(self):
		self.clinvar_path = "/groups/clinical/projects/editability/clinvar/"
		self.CRISPRpath = "/groups/clinical/projects/editability/CRISPR_hg38/"
		self.HGVSlookup_path = "/groups/clinical/projects/editability/clinvar/HGVSlookup.csv"
		self.HGVSlookup_table = None
		self.guide_tab = None


	def updateTables(self):
		'''
		checks date on clinvar table and if > month than runs aggClinvar.py
		updateClinVarHPA(clinvar_ftp = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
                 out_dir = "/groups/clinical/projects/editability/clinvar/",
                 HPA_file="/groups/clinical/projects/editability/HPA_RNA_TissueExp/proteinatlas.tsv")
		'''
		pass

	def get_HGVStable(self):
		'''
		Get HGVS Lookup table - a table to find the chrom a HGVS ID resides
		'''
		hgvs_tab = pd.read_csv(self.HGVSlookup_path)
		self.HGVSlookup_table = hgvs_tab
		return hgvs_tab


	def get_Guidetable(self,chrom,pos,win_size):
		bbfile = f"{self.CRISPRpath}chr{str(chrom)}.bb"
		f = bbi.open(bbfile)
		if pos == None:
			end = list(bbi.chromsizes(bbfile).values())
			df = f.fetch_intervals(chrom=f"chr{chrom}", start=1, end=end)
		else:
			win_dif = win_size[1]-win_size[0]
			start1,end1 = pos-win_size[0], pos-win_size[0] +23
			start2,end2 = (pos-win_size[0]) + win_dif, pos-win_size[0] +23+ win_dif
			df1 = f.fetch_intervals(chrom=f"chr{chrom}", start=start1, end=start1)
			df2 = f.fetch_intervals(chrom=f"chr{chrom}", start=start2, end=start2)
			df = pd.concat([df1, df2], ignore_index=True).drop_duplicates()

		self.guide_tab = df
		return df

	def get_ClinVartable(self,chrom):
		df = pd.read_csv(f"{self.clinvar_path}{chrom}_variant.txt")
		return df

	def getChrom(self,hgvs_id):
		hgvs_tab = self.get_HGVStable()
		hgvs_list = hgvs_tab["HGVS"]
		try:
			chrom = hgvs_tab['Chr'].iloc[list(hgvs_list).index(hgvs_id.split(".")[0])]
			return chrom
		except IndexError:
			raise "No HGVSs were found for this search"


#####
# Data Query Functions
#####
def aggData(guides,vardf):
	n_guides = guides.shape[0]
	print("------------------------------------------------")
	print("HGVS ID: ",vardf.Name.iloc[0])
	print("SNV Type: ", vardf.Type.iloc[0])
	print("Gene Name: ", vardf.GeneSymbol.iloc[0])
	print("  >>>>>> Phenotype <<<<<<   ")
	print("Listed Phenotypes: ", vardf.PhenotypeList.iloc[0])
	print("Tissue Enrichment: ", vardf['RNA tissue cell type enrichment'].iloc[0])
	print("Molecular Function: ", vardf['Molecular function'].iloc[0])
	print("  >>>>>> Guides Found <<<<<< ")
	cnt = 1
	variant_pos = vardf["PositionVCF"].iloc[0]
	'''
		
		CCC-TTAGCTTCGCCGACAACCCC
		GGG-AATCGAAGCGGCTGTTGGGG
	'''
	for i in range(n_guides):
		print(f"  ---- Guide{cnt} -----")
		pos = variant_pos - guides.start.iloc[i]
		print("        ","".join(" " if p != pos else "*" for p in range(23)))
		print("guideSeq: ", guides.field8.iloc[i])
		print("PAM: ", guides.field9.iloc[i])
		x = guides.field10.iloc[i]
		scores = re.findall("(\w+:[()%0-9]+)",x.replace(" ",""))
		for s in scores:
			print(f"  {s}")
		cnt+=1

def query_HGVS(hgvs_id,win_size=[4,8]):
	# test   hgvs_id = 'NM_000152.5:c.271G>A'
    	dh = DataHandler()

	if hgvs_id.startswith('NM'):

		try:
			prefix = re.search("NM_\d+.",hgvs_id).captures()[0]
		except:
			print("Be sure hgvs_id contains the prefix ex: NM_0000")

		chrom = dh.getChrom(hgvs_id)
		vardf = dh.get_ClinVartable(chrom)
		vardf = vardf.loc[vardf["Name"].str.startswith(prefix)]
		try:
			var_desc = re.search("\d+\w+>\w+", hgvs_id).captures()[0]
			try:
				vardf = vardf.loc[vardf["Name"].str.contains(var_desc)]
			except IndexError:
				raise "No HGVSs were found for this search"
		except:
			print("Be sure hgvs_id contains the variant description ex: 271G>A")

		variant_pos = vardf["PositionVCF"].iloc[0]

		guides = dh.get_Guidetable(chrom,variant_pos,win_size)
				aggData(guides, vardf)


def queryGuides(query,qtype):
	if qtype == 'HGVS':
		query_HGVS(hgvs_id=query, win_size=[4, 8])

queryGuides(query,qtype)


		return guides, vardf



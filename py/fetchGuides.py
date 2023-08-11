import pandas as pd
import bbi
import regex as re

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
			start = pos-win_size[0]
			end = pos+23+(win_size[1]-win_size[0])
			df = f.fetch_intervals(chrom=f"chr{chrom}", start=start, end=start)
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
				vardf = vardf.loc[df["Name"].str.contains(var_desc)]
			except IndexError:
				raise "No HGVSs were found for this search"
		except:
			print("Be sure hgvs_id contains the variant description ex: 271G>A")

		variant_pos = vardf["PositionVCF"].iloc[0]

		guides = dh.get_Guidetable(chrom,variant_pos,win_size)


		return guides, vardf



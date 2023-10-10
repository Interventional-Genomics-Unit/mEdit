# Native Modules
# import regex as re
# import sys
import os
from datetime import date  # datetime, date
# Installed Modules
import pandas as pd
# from Bio.Seq import Seq
# from Bio import motifs
# Project Modules
from dataH import DataHandler, get_seqinfo

###############
# Main Script with Fetch_Guides Class for running pipeline
###############

# from py.validate import Validator


def set_export(outdir):
	# Create outdir if inexistent
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	return outdir


class Fetch_Guides:

	def __init__(self,
				 queries: list,
				 qtype: str,
				 editor: str,
				 BEmode: str,
				 resultsfolder: str,
				 datadir: str,
				 fasta_path: str,
				 **kwargs):
		"""
		:param queries: list of query terms
		:param qtype: 'hgvs' or 'coord' #coord is not function yet
		:param editor: 'all', 'custom'(custom must contain kwargs, see below),
		or selected list or str() of the editor choices
		:param BEmode: 'off','default','all', or select BE editor for base editor choices below
		:param kwargs:

		** if 'custom' selcted as editor in kwargs must include pam, pamISFirst, window_size (optional:name)
		** if not custom or
		"""

		# -----------------User Inputs-------------------- #
		self.queries = queries
		self.qtype = qtype
		self.editor = editor
		self.BEmode = BEmode
		self.kwargs = kwargs

		# input paths
		self.datadir = datadir
		self.processed_tables = f"{self.datadir}/processed_tables"  # folder with cleaned clinvar/hpa tabs
		self.HGVSlookup_path = f"{self.processed_tables}/HGVSlookup.csv"
		self.fasta_path = fasta_path

		# output folder
		self.resultsfolder = resultsfolder

		# ---------------libraries and keys--------------------#
		# [editor]: pam, pamISfirst, win_size, guidelen, scoring, notes/altnames
		self.editor_choices = ['spCas9', 'saCas9', 'spG', 'SpRY-HighE', 'LbCpf1',
							   'scCas9', 'stCas9', 'iSpyMacCas9', 'CasX', 'Cpf1']

		self.BE_choices = ['BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-Gam', 'YE1-BE3', 'EE-BE3', 'YE2-BE3',
						   'YEE-BE3', 'VQR-BE3', 'VRER-BE3', 'SaBE3', 'SaBE4', 'SaBE4-Gam', 'Sa(KKH)-BE3', 'xBE3',
						   'eA3A-BE3',
						   'A3A-BE3', 'BE-PLUS', 'ABE7.9', 'ABE7.10', 'xABE,ABESa', 'VQR-ABE', 'VRER-ABE',
						   'Sa(KKH)-ABE']

		self.editor_pamlib = {'spCas9':('NGG', False,20,[10,24], 'Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
                              'saCas9': ('NNGRRT', False, 21,[10,24], 'Cas9 S. Aureus 21 base guide'),
                              'spG':('NGN',False,20,[10,24],'20bp-NGN - SpG'),
                              'SpRY-HighE': ('NRN',False,20,[10,24],'High Efficiency Pam'),
                              'scCas9':('NNGT',False,20,[10,24],'20bp-NNGT - Cas9 S. canis - high efficiency PAM, recommended'),
                              'stCas9':('NNAGAA',False,20,[10,24],'Cas9 S. Thermophilus'),
                              'iSpyMacCas9':('NAA',False,20,[10,24],''),
                              'CasX':('TTCN',True,20,[0,7],'Cas12e'),
                              'AsCas12a':('TTTV',True,23,[0,7],'TTT(A/C/G)-23bp - Cas12a (Cpf1)'),
                              'LbCas12a': ('TTTV', True, 23,[0,7], 'LbCpf1'),
                              'Cas12c1': ('TG', True, 23, [0, 7], 'C2c3'),
							  }

		## ------------Defaults and settings------------------##
		##configure editor options
		self.search_params = self.configure_search_params()

		# configure BE options
		if self.BEmode != 'off':
			self.BE_search_params = self.set_BE_params()

		# ---------------Ouputs--------------------------##
		self.all_variant = pd.DataFrame()
		self.all_gene = pd.DataFrame()
		self.all_guides = {}
		self.all_BE = {}

		self.found_genes = []

	def set_guidelen(self, guidelen):
		self.guidelen = guidelen

	def set_win_size(self, win_size):
		self.win_size = win_size

	def configure_search_params(self):
		"""
		set paramteres for the selected editor or editors(not BE editors)
		"""
		# search for all guides
		if 'all' == self.editor:
			self.search_params = self.editor_pamlib

		# search for selected subset
		if type(self.editor) is list:
			self.search_params = {}
			for e in self.editor:
				self.search_params[e] = self.editor_pamlib[e]

		# else use single set parameters
		else:
			# default - spCas9 params
			self.win_size = [4, 8]
			self.guidelen = 20
			self.scoring = 'doench'
			self.pamISfirst = False
			self.pam = 'NGG'
			self.name = 'spCas9'
			self.notes = 'none'

			# set custom editor params
			if 'custom' == self.editor:
				self.search_params = self.set_params(self.kwargs)

			# select a single editor
			if self.editor in self.editor_choices:
				self.search_params = {self.editor: self.editor_pamlib[self.editor]}

		return self.search_params

	def set_params(self, kwargs):
		# opts = ['pam', 'pamISfirst', 'guidelen','win_size', 'name','scoring']

		if 'pam' in kwargs.keys():
			self.pam = kwargs['pam']
		if 'pamISfirst' in kwargs.keys():
			self.pamISfirst = kwargs['pamISfirst']
		if 'guidelen' in kwargs.keys():
			self.guidelen = kwargs['guidelen']
		if 'name' in kwargs.keys():
			self.name = kwargs['name']
		if 'win_size' in kwargs.keys():
			self.win_size = kwargs['win_size']
		if 'notes' in kwargs.keys():
			self.notes = kwargs['notes']

		params = {self.name: (self.pam, self.pamISfirst, self.win_size, self.guidelen)}

		return params

	def set_multi_params(self):
		# TODO: make automatic function using configureMultipams
		pass

	def set_BE_params(self):
		# sets base editor search params, each key is a list of 2 or more; refernce seq search params,
		# then any set that follows starts with the conversion (ex. 'AG' is A --> G) and then the base editors that have the same params

		BE_lib = {'spCas9-def': [('NGG', False, 20, [4, 8]), ('CT', 'BE3', 'BE4', 'BE4max', 'BE4-Gam'),
								 ('AG', 'ABE7.9', 'ABE7.10', 'ABEmax')],
				  'spCas9-YE': [('NGG', False, 20, [5, 6]), ('CT', 'YE1-BE3', 'EE-BE3', 'YE2-BE3', 'YEE-BE3')],
				  'spCas9-VQR': [('NGA', False, 20, [4, 8]), ('CT', 'VQR-BE3'), ('AG', 'VQR-ABE')],
				  'spCas9-VRER': [('NGCG', False, 20, [4, 8]), ('CT', 'VRER-BE3'), ('AG', 'VRER-ABE')],
				  'saCas9-KKH': [('NNGRRT', False, 20, [3, 12]), ('CT', 'SaBE3', 'SaBE4', 'SaBE4-Gam,Sa(KKH)-BE3'),
								 ('AG', 'ABESa', 'Sa(KKH)-ABE')]}

		if self.BEmode == 'default':
			self.BE_search_params = {'spCas9-def': BE_lib['spCas9-def']}

		if self.BEmode == 'all':
			self.BE_search_params = BE_lib

		if self.BEmode in self.BE_choices:
			for k, v in BE_lib.values():
				if self.BEmode in v[1]:
					self.BE_search_params = [BE_lib[k]('CT', self.BEmode)('AT')]
				if self.BEmode in v[2]:
					self.BE_search_params = [BE_lib[k]('CT')('AT', self.BEmode)]

		return self.BE_search_params

	def write_guide_csv(self, guides, gtype):
		df = pd.DataFrame(guides)
		nameout = 'BaseEditors_found.csv' if gtype == 'BE' else 'Guides_found.csv'
		datenow = date.today().strftime('%Y-%m-%d')
		out = f'{datenow}_{nameout}.csv'
		df.to_csv(f"{self.resultsfolder}/{out}", index=False)
		return df

	def add_clininfo(self):
		clininfo = pd.DataFrame()
		ids = set(self.all_guides['HGVS_ID'])
		chroms = self.all_guides['Chr']
		for ch in set(chroms):
			tempdf = pd.read_csv(f"{self.processed_tables}/{ch}_variant.txt")
			tempdf = tempdf.loc[tempdf['HGVS_ID'].isin(list(ids))]
			clininfo = pd.concat([clininfo, tempdf])
		self.all_variant = clininfo[
			['HGVS_ID', 'GeneSymbol', 'Chr', 'PositionVCF', 'Strand', 'RefAlleleVCF', 'AltAlleleVCF',
			 'AlleleID', 'Type', 'GeneID', 'HGNC_ID', 'ClinicalSign', 'ClinSigSimple', 'LastEval',
			 'RS#(dbSNP)', 'nsv/esv (dbVar)', 'RCVaccession', 'PhenoList', 'Origin', 'OriginSimple',
			 'ChrAccession', 'Cytogenetic', 'ReviewStatus', 'NumberSubmitters', 'Guidelines', 'TestedInGTR',
			 'SubmitterCategories', 'VariationID', 'OMIM', 'IDs']]
		if 'Protein class' in clininfo.columns:
			self.all_gene = clininfo[['GeneSymbol', 'Chr', 'Start', 'End', 'ChrID',
			                          'TranscriptID', 'ProteinID', 'Ensembl_Gene', 'Protein class', 'Gene description',
			                          'Biological process', 'Molecular function', 'Uniprot',
			                          'Disease involvement', 'RNA tissue specificity', 'RNA tissue specific nTPM',
			                          'RNA tissue distribution', 'RNA tissue cell type enrichment']]
		else:
			self.all_gene = clininfo[['GeneSymbol', 'Chr', 'Start', 'End', 'ChrID',
			                          'TranscriptID', 'ProteinID', 'Ensembl_Gene']]
		datenow = date.today().strftime('%Y-%m-%d')

		self.all_gene.to_csv(f'{self.resultsfolder}/{datenow}_Gene_Report.csv',index = False)
		self.all_variant.to_csv(f'{self.resultsfolder}/{datenow}_Variant_Report.csv',index = False)

	def run_FetchGuides(self):
		# information needed for guide aquisition
		variantseq_dict = get_seqinfo(self.queries, self.qtype, self.datadir)
		for hgvs_id, data in variantseq_dict.items():
			dh = DataHandler(hgvs_id, data, self.fasta_path)
			if self.BEmode != 'off':
				guides, BEguides = dh.get_Guides(self.search_params, self.BE_search_params)
			else:
				guides, BEguides = dh.get_Guides(self.search_params)
			if len(BEguides['gRNA']) > 0:
				if len(self.all_BE.keys()) == 0:
					for k, v in BEguides.items():
						self.all_BE[k] = v
				else:
					for k, v in BEguides.items():
						self.all_BE[k] += v
			if len(guides['gRNA']) > 0:
				if len(self.all_guides.keys()) == 0:
					for k, v in guides.items():
						self.all_guides[k] = v
				else:
					for k, v in guides.items():
						self.all_guides[k] += v
				print(len((guides['gRNA'])), ' guides found for ', hgvs_id)
			else:
				print(f"No guides found for the query {hgvs_id}")
		guidedf = self.write_guide_csv(self.all_guides, gtype='guides')
		BEdf = self.write_guide_csv(self.all_BE, gtype='BE')
		self.add_clininfo()
		return self.all_variant, self.all_gene, guidedf, BEdf


def main():
	# SNAKEMAKE IMPORTS
	#   Inputs
	input_file = str(snakemake.input.query_manifest)
	fasta_path = str(snakemake.input.assembly_path)
	#   Outputs
	resultsfolder = set_export(str(snakemake.output))
	#   Params
	datadir = str(snakemake.params.support_tables)
	# Paths---------------------------
	# input_file = "/groups/clinical/projects/editability/medit_queries/medit_test/test_in/hgvs_test_queries.csv"
	# resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"
	# datadir = "/groups/clinical/projects/editability/tables/"
	# fasta_path = "/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa"

	# Input Extraction-------------------
	df = pd.read_csv(input_file)
	queries = list(df.iloc[:, 0])
	qtype = 'hgvs'
	BEmode = 'default'
	editor = 'all'

	# queries += ['NM_000152.5(GAA):c.271G>T']

	# Report processed input variables
	print(f"""
	Currently running fetchGuides.py
	INPUT VARIABLES:
		Queries:\n{queries}
		Query Type: {qtype}
		BEmode: {BEmode}
		editor: {editor}
	PATH TO REFERENCE:
		-> {fasta_path}
	SUPPORT DATA DIRECTORY:
		-> {datadir}
	""")
	# Get query items
	fg = Fetch_Guides(queries, qtype, editor, BEmode, resultsfolder, datadir, fasta_path)
	all_clin_info, all_gene, all_guides, all_BE = fg.run_FetchGuides()


if __name__ == "__main__":
	main()

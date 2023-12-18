# Native Modules
import gzip
import zlib

# import regex as re
import os
import re
from zlib import error
# Installed Modules
import pandas as pd
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
import pickle
# Project Modules
from dataH import DataHandler
from annotate import get_refseq_entry, find_snvseq_info, find_transcript_info

###############
# Main Script with Fetch_Guides Class for running pipeline
###############


def set_export(outdir):
	# Create outdir if inexistent
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	return outdir


class Fetch_Guides:

	def __init__(self,
	             queries: list,
	             qtype: str,
	             editor_request: str | list,
	             be_request: str | list,
	             editors: dict,
	             base_editors: dict,
	             datadir: str,
	             fasta_path: str,
	             annote_path: str,
	             **kwargs):
		"""
		:param queries: list of query terms, either in hgvs format - 'NM_000518.5:c.114G>A' or coords 'chr11:5226778C>T' (COORDS ALLELES MUST BE PLUS STRAND!!)
		:param qtype: 'hgvs' or 'coord'
		  --> if 'hgvs', providing the coordinates in the kwargs with 'hgvscoord' can reduce processing time
		  --> hgvs assumes the query is already in clinvar and will generate a variant report with the gene report,
		  --> if 'coord' then just gene report is created
		:param editor_request: 'clinical', 'custom', name (list/str from editor choices)
		--> custom must contain kwargs - pam, pamISFirst,guidelen (optional:name,win_size)
		:param be_request: 'off','default','all', or select BE editor for base editor choices below
		:param editors: Dictionary containing information on the current set of editors supported by mEdit
		:param base_editors: Dictionary containing information on the current set of base editors supported by mEdit
		:param genome: genome used
		:param datadir: folder where tables and pre-computed data live
		:param fasta_path: *Unsure using chromsome seperate files right now but unsure if this will be permenant
		:param kwargs: 'hgvscoord' , 'clin_report','gene_report'

		"""
		##-----------------User Inputs--------------------##
		if qtype == 'hgvs':
			self.queries = self.validate_hgvs(queries)
		if qtype == 'coord':
			self.queries = self.validate_coord(queries)
		self.qtype = qtype
		self.editor_request = editor_request
		self.be_request = be_request
		self.hgvscoord = None
		self.clin_report = True
		self.gene_report = True
		self.kwargs = kwargs
		self.editor_lib = editors
		self.be_lib = base_editors

		if 'hgvscoord' in kwargs.keys():
			self.hgvscoord = self.validate_coord(kwargs['hgvscoord'])  # 'chr11:5226778C>T'
		if 'gene_report' in kwargs.keys():
			self.gene_report = kwargs['gene_report']
		if 'clin_report' in kwargs.keys():
			self.clin_report = kwargs['clin_report']

		# input paths/folders
		self.processed_tables = f"{datadir}/processed_tables"  # folder with cleaned clinvar/hpa tabs
		self.HGVSlookup_path = f"{self.processed_tables}/HGVSlookup.csv"
		self.fasta_path = fasta_path
		self.annote_path = annote_path

		# other variables
		self.snv_info = {}  # {chrom: (id,snv_pos,ref,alt)}

		##---------------libraries and keys--------------------##
		self.editor_choices = list(editors.keys())
		self.BE_choices = list(base_editors.keys())

		## ------------Defaults and settings------------------##
		##configure editor options
		self.search_params = self.configure_search_params()

		# configure BE options
		if self.be_request != 'off':
			self.BE_search_params = self.set_BE_params()
		# ---------------Flags--------------------------#
		self.clininfo_flag = False

		# ---------------Ouputs--------------------------#
		self.all_variant = pd.DataFrame()
		self.all_gene = pd.DataFrame()
		self.all_guides = {}
		self.all_BE = {}


	def configure_search_params(self):
		"""
		set paramteres for the selected editor or editors(not BE editors)
		"""
		# search for all guides
		if self.editor_request == 'clinical':
			search_params = self.editor_lib['clinical']
		# set custom editor params
		elif self.editor_request == 'custom':
			search_params = self.set_params(self.kwargs)

		# search for selected subset
		# TODO: The guide_prediction.py needs to ingest this correctly
		elif type(self.editor_request) is list:
			search_params = {}
			for e in self.editor_request:
				search_params[e] = self.editor_lib[e]

		# else use single set parameters
		else:
			if self.editor_request in self.editor_choices:
				search_params = {self.editor_request: self.editor_lib[self.editor_request]}

			else:
				print('Please choose a name(s) found in the editor name choices')

		print(f'Editor(s) set to: {[x for x in search_params.keys()]}')
		return search_params

	def set_params(self, kwargs):
		name = 'custom'

		try:
			pam = kwargs['pam']
			pamISfirst = kwargs['pamISfirst']
			guidelen = kwargs['guidelen']
			if pamISfirst == False:
				win_size = -2
			else:
				win_size = guidelen - 2
		except KeyError:
			print("custom editor selection MUST include a minimum of kwargs = pam, pamISFirst,guidelen")

		if 'name' in kwargs.keys():
			name = kwargs['name']
		if 'win_size' in kwargs.keys():
			win_size = kwargs['win_size']

		params = {name: (pam, pamISfirst, win_size, guidelen, '')}

		return params

	def set_BE_params(self):
		# sets base editor search params, each key is a list of 2 or more; refernce seq search params,
		# then any set that follows starts with the conversion (ex. 'AG' is A --> G) and then the base editors that have the same params

		if self.be_request == 'default':
			self.BE_search_params = self.be_lib['default']
		elif self.be_request == 'all':
			self.BE_search_params = self.be_lib['all']

		# TODO: Setup BE user-defined list processing
		elif type(self.be_request) is list:
			self.BE_search_params = {}
			for e in self.be_request:
				self.BE_search_params[e] = self.be_lib['all'][e]

		else:
			if self.be_request not in self.BE_choices:
				print('That is not a valid Base Editor')
				print(f'please choose from {self.BE_choices}')
			else:
				for k, v in self.be_lib['all'].values():
					if self.be_request in v[1][-1]:
						self.BE_search_params = {self.be_request: self.be_lib['all'][k][0:2]}

					if len(v) == 3:
						if self.be_request in v[2][-1]:
							self.BE_search_params = {self.be_request: self.be_lib['all'][k][0] + self.be_lib['all'][k][2]}

		return self.BE_search_params

	def write_gsearch_params(self, outfile):
		# writes pickle of selected guide search params for later use in process_genome
		# 'editor', 'pam', '5prime_pam','guide_length', 'DSB site', 'notes'
		with open(outfile, 'ab') as gfile:
			pickle.dump(self.search_params, gfile)


	def write_snv_site_info(self, outfile):
		'''
		#writes pickle of SNV site info for later use in process genome
		#query, tid, eid, strand, ref, alt, feature_annotation, extracted_seq, codons, coord
		'''
		with open(outfile, 'ab') as sfile:
			pickle.dump(self.snv_info, sfile)

	def write_guide_csv(self, guides, outfile):
		df = pd.DataFrame(guides)
		if 'On-Target Efficiency Score' in df.columns:
			temp = df[df['On-Target Efficiency Score'] != '-'].sort_values(by='On-Target Efficiency Score', ascending=False)
			df = pd.concat([temp, df[df['On-Target Efficiency Score'] == '-']]).reset_index(drop=True)
		df['Guide_ID'] = [y + str(x) for x, y in zip(list(df.index), list(df['Guide_ID']))]
		df.to_csv(outfile, index=False)
		return df

	def add_clininfo(self, gene_out, variant_out):
		if not self.clininfo_flag:
			# TAYLOR: Here we can probably provide something more informative.
			#   Keeping it as a placeholder
			self.all_gene.to_csv(gene_out)
			print("GENES AND VARIANT TABLES ARE UNAVAILABLE")
			return
		all_tids = []
		for ch, data in self.snv_info.items():
			all_tids += [d[1] for d in data]

			if self.qtype == 'hgvs':
				tempvar = pd.read_csv(f"{self.processed_tables}/variant_tables/{ch}_variant.txt")
				tempvar = tempvar.loc[tempvar['HGVS_Simple'].isin(list(self.queries))]
				self.all_variant = pd.concat([self.all_variant, tempvar])
		#TODO: The path to gene_tables will ideally be imported as a snakemake variable
		tempgene = pd.read_csv(f"{self.processed_tables}/gene_tables/gene_tables.csv")
		self.all_gene = tempgene.loc[tempgene['TranscriptID'].isin(list(all_tids))]

		# datenow = date.today().strftime('%Y-%m-%d')

		# gene_out = f"{self.resultsfolder}/Gene_Report.csv"
		print(f"\n READY TO PRINT GENE OUT TO: {gene_out}\n ")
		self.all_gene.to_csv(gene_out, index=False)

		if self.qtype == 'hgvs':
			# variant_out = f"{self.resultsfolder}/Variant_Report.csv"
			print(f"\nREADY TO PRINT VARIANT OUT TO: {variant_out}\n")
			self.all_variant.to_csv(variant_out, index=False)

	@staticmethod
	def extract_seqs(searchseq, pos, alt, window):
		"""
		extracts the sequence +/-30bp surrounding a SNV then swaps ref for alt allele
		"""
		extracted_seq = str(searchseq[pos - window:pos + window])
		extracted_seq = Seq(extracted_seq[0:window] + alt + extracted_seq[window + 1:]).upper()
		return extracted_seq

	def get_chroms(self,queries):
		hgvs_tab = pd.read_csv(self.HGVSlookup_path)
		q_prefixes = [x.split(':')[0] for x in queries]
		chroms = set(hgvs_tab.loc[hgvs_tab['TranscriptID'].isin(q_prefixes), 'Chr'])
		return chroms

	def fetch_query_info(self):
		# Gets Transcript info
		global term
		snv_info = {}
		window = 50

		# If quering by HGVSID with no other info then need to get chromsome/location/alt/ref
		if self.qtype == 'hgvs' and self.hgvscoord is None:
			print("Looking up HGVS in Clinvar.......")
			chroms = self.get_chroms(self.queries)

			for ch in chroms:
				df = pd.read_csv(f"{self.processed_tables}/variant_tables/{ch}_variant.txt", low_memory=False)
				gadf = df.loc[df['HGVS_Simple'].isin(self.queries)]
				snv_info[ch] = gadf[['HGVS_Simple', 'PositionVCF', 'RefAlleleVCF', 'AltAlleleVCF']].to_dict('tight')[
					'data']

		# Else All information is given to find transcript info
		else:
			coords = self.queries if self.qtype == 'coord' else self.hgvscoord

			coord_fmt = r'chr[0-9MTXY]*:(\d*)([ATCG]{1})\>([ATCG]{1})'

			print(f"\n BUG INSPECTION \n Coords: {coords}\n Queries: {self.queries}\n")
			print(f"PREMISSAS:\n Qtype: {self.qtype}\n Hgvs Coord: {self.hgvscoord}")
			for x in range(len(self.queries)):
				print(f" Current Query: {x}\n")
				ch = coords[x].split(':')[0].replace('chr', '')
				if ch not in snv_info.keys():
					snv_info[ch] = []
				snvpos, alt, ref = list(re.search(coord_fmt, coords[x]).groups())
				snv_info[ch].append([self.queries[x], int(snvpos), alt, ref])

		self.snv_info = snv_info

		print("Gathering Variant Genomic Info.......")

		for ch, data in snv_info.items():  # find transcript info
			# === Transitioning SeqIO.read to direct import of Pickled SeqRecord Objects ===
			# chr_fasta_path = self.fasta_path.replace('.fa.gz', f'_chr{str(ch)}.fa.gz')
			chr_fasta_path = f"{self.fasta_path}/chr{str(ch)}.pkl"
			try:
				print(f"Finding transcripts information: Assessing {chr_fasta_path}")
				# fasta = SeqIO.read(gzip.open(chr_fasta_path, 'rt'), 'fasta')
				with open(chr_fasta_path, 'rb') as pfile:
					fasta = pickle.load(pfile)
				new_data = []
			except FileNotFoundError:
				print(f"The file {chr_fasta_path} was not found. Please regenerate background data and check the target directory")
				continue
			except pickle.UnpicklingError:
				print(f"The file {chr_fasta_path} is not in the correct format. Please regenerate background data")
				continue

			for d in data:
				query, snvpos, ref, alt = d
				print(query, ref, alt, snvpos)
				if self.qtype == 'hgvs':  # pull refseqID from HGVS and search transcript by this
					term = query.split(':')[0]
				if self.qtype == 'coord':  # else use coordsinates to search trancript
					term = f"chr{str(ch)}:{str(snvpos)}-{str(snvpos)}"
				entry, tid_info = find_transcript_info(term=term, fasta=fasta, annote_path = self.annote_path)

				if entry is not None:
					exons, tx_seq, cds = tid_info
					t_snvpos = int(snvpos) - int(entry['txStart'])
					extracted_seq = self.extract_seqs(searchseq=tx_seq, pos=t_snvpos - 1, alt=alt, window=window)
					if len(extracted_seq) != window * 2:
						# if flanking or utr the extracted seq needs to come from the chromosome file
						extracted_seq = self.extract_seqs(searchseq=fasta.seq[snvpos - 100:snvpos + 100], pos=t_snvpos - 1,
												alt=alt,
												window=window)

					feature_annotation, codons = find_snvseq_info(extracted_seq,snvpos,exons, cds, entry,window)
					strand = entry['strand']
					tid, eid = entry['tid'], entry['eid']

				else:
					feature_annotation = 'undetermined/non-coding'
					codons = 'None'
					strand = '+'
					extracted_seq = self.extract_seqs(searchseq=fasta.seq, pos=snvpos, alt=alt, window=window)
					tid, eid = term, '-'
				new_data.append(
					[query, tid, eid, strand, ref, alt, feature_annotation, extracted_seq, codons,
					 f"chr{str(ch)}:{str(snvpos)}"])

				print('Query term & annoation:', query, feature_annotation)
			snv_info[ch] = new_data

		self.snv_info = snv_info

	@staticmethod
	def validate_hgvs(queries):
		'''
		standardizes input hgvs and checks formating
		'''
		rprefix = r"((N(M|G|C|R)_[\d.]*)|(m))"
		rsuffix = r"(:(c|m|g|n)\S*)"
		validated_queries = []
		for q in set(queries):
			if re.search(rsuffix, q) and re.search(rprefix, q):
				validated_queries.append(re.search(rprefix, q).groups()[0] + re.search(rsuffix, q).groups()[0])
			else:
				print(q)
		n = len(validated_queries)
		print(f'{n} out of {len(queries)} HGVS IDs validated')
		if n == 0:
			print('Query are not in the correct HGVS Format')
		return validated_queries

	@staticmethod
	def validate_coord(queries):
		'''
		standardizes input coordinate and checks formatting
		'''
		# q = 'chr11:5226778C>T'
		coord_fmt = r'(chr[0-9]*:\d*(A|T|C|G)>(A|T|C|G))'
		validated_queries = []
		for q in set(queries):
			if re.search(coord_fmt, q):
				validated_queries.append(re.search(coord_fmt, q).groups()[0])
			else:
				print(q)
		n = len(validated_queries)
		print(f'{n} out of {len(queries)} Coordinates IDs validated')
		if n == 0:
			print('Query are not in the correct Coordinate + allele Format')
		return validated_queries

	def run_FetchGuides(self, outfile_guides, outfile_be_guides):
		global dh, query
		self.fetch_query_info()
		print('Finding Guides.....')
		for ch, data in self.snv_info.items():

			for d in data:
				try:
					query, tid, eid, strand, ref, alt, feature_annotation, extracted_seq, codons, coord = d
				except ValueError:
					print(f"WARNING: The query below has the wrong number of values to unpack. Needs further investigation:\n{d}")
					continue
				dh = DataHandler(query, strand, ref, alt, feature_annotation, extracted_seq, codons, coord)

				if self.be_request != 'off':
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

						print(len((guides['gRNA'])), ' guides found for ', query)
				else:
					print(f"No guides found for the query {query}")

		guidedf, BEdf = None, None

		if len(self.all_guides.keys()) != 0:
			guidedf = self.write_guide_csv(self.all_guides, outfile_guides)
			self.clininfo_flag = True
			# self.add_clininfo()

		if len(self.all_BE.keys()) != 0:
			BEdf = self.write_guide_csv(self.all_BE, outfile_be_guides)
			self.clininfo_flag = True

		return {'all_variant': self.all_variant,
		        'all_gene': self.all_gene,
		        'guide_table': guidedf,
		        'BE_table': BEdf}


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	input_file = str(snakemake.input.query_manifest)
	fasta_path = str(snakemake.input.assembly_path)
	# === Outputs ===
	guides_report = str(snakemake.output.guides_report_out)
	guide_search_params_path = str(snakemake.output.guide_search_params)
	snv_site_info_path = str(snakemake.output.snv_site_info)
	# === Params ===
	resultsfolder = set_export(str(snakemake.params.main_out))
	gene_report = f"{resultsfolder}/{str(snakemake.params.gene_report)}"
	variant_report = f"{resultsfolder}/{str(snakemake.params.variant_report)}"
	be_report = f"{resultsfolder}/{str(snakemake.params.be_report)}"
	#   == Processed tables branch ==
	datadir = str(snakemake.params.support_tables)
	annote_path = str(snakemake.params.annote_path)
	# == Editor Parameters
	editors_path = str(snakemake.params.editors)
	base_editors_path = str(snakemake.params.base_editors)
	#   == Run Parameters ==
	qtype = str(snakemake.params.qtype)
	be_request = str(snakemake.params.be_request)
	editor_request = str(snakemake.params.editor_request)
	# == DEBUG BLOCK ==
	# qtype = 'hgvs'
	# BEmode = 'off'
	# editor = 'all'
	# == == ==

	# == Input Setup ==
	df = pd.read_csv(input_file)
	queries = list(df.iloc[:, 0])

	# == Editors / BEs Setup ==
	with open(editors_path, 'rb') as edfile:
		editors = pickle.load(edfile)
	with open(base_editors_path, 'rb') as befile:
		base_editors = pickle.load(befile)

	# == Report processed input variables ==
	print(f"""
	Currently running fetchGuides.py
	INPUT VARIABLES:
		Queries:\n{queries}
		Query Type: {qtype}
		be_request: {be_request}
		editor_request: {editor_request}
	PATH TO REFERENCE:
		-> {fasta_path}
	SUPPORT DATA DIRECTORY:
		-> {datadir}
	OUTPUTS TO:
		--> {resultsfolder}
	""")

	# == Get query items ==
	fg = Fetch_Guides(queries,
	                  qtype,
	                  editor_request,
	                  be_request,
	                  editors,
	                  base_editors,
	                  datadir,
	                  fasta_path,
	                  annote_path
	                  )
	# == Set up object and run core methods ==
	exports = fg.run_FetchGuides(guides_report, be_report)

	# == Export Intermediate files ==
	fg.write_snv_site_info(snv_site_info_path)
	fg.write_gsearch_params(guide_search_params_path)

	# == Export Variant and Gene tables ==
	fg.add_clininfo(gene_report, variant_report)


if __name__ == "__main__":
	main()

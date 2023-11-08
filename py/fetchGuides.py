# Native Modules
import gzip
import time
import regex as re
import os
from datetime import date
# Installed Modules
import pandas as pd
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
import pickle
# Project Modules
from dataH import DataHandler


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
	             editor: str | list,
	             BEmode: str | list,
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
		:param editor: 'all', 'custom', selected list or str() of the editor choices
		--> custom must contain kwargs, see below
		:param BEmode: 'off','default','all', or select BE editor for base editor choices below
		:param genome: genome used
		:param datadir: folder where tables and pre-computed data live
		:param fasta_path: *Unsure using chromsome seperate files right now but unsure if this will be permenant
		:param kwargs: 'hgvscoord' , 'clin_report','gene_report'

		** if 'custom' selcted as editor in kwargs must include pam, pamISFirst, window_size (optional:name)
		"""
		##-----------------User Inputs--------------------##
		self.notes = None
		self.name = None
		self.pam = None
		self.pamISfirst = None
		self.scoring = None
		self.win_size = None
		self.guidelen = None
		if qtype == 'hgvs':
			self.queries = self.validate_hgvs(queries)
		if qtype == 'coord':
			self.queries = self.validate_coord(queries)
		self.qtype = qtype
		self.editor = editor
		self.BEmode = BEmode
		self.hgvscoord = None
		self.clin_report = True
		self.gene_report = True
		self.kwargs = kwargs

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
		# [editor]: pam, pamISfirst, win_size, guidelen, scoring, notes/altnames
		self.editor_choices = ['spCas9', 'saCas9', 'spG', 'SpRY-HighE', 'scCas9',
		                       'stCas9', 'iSpyMacCas9', 'CasX', 'AsCas12a', 'LbCas12a', 'Cas12c1']

		self.BE_choices = ['BE1', 'BE3', 'BE2', 'HF-BE3', 'BE4', 'BE4max', 'BE4-Gam', 'YE1-BE3', 'EE-BE3', 'YE2-BE3',
		                   'YEE-BE3', 'VQR-BE3', 'VRER-BE3', 'SaBE3', 'SaBE4', 'SaBE4-Gam', 'Sa(KKH)-BE3', 'xBE3',
		                   'Target-AID',
		                   'ABE7.9', 'ABE7.10', 'xABE,ABESa', 'VQR-ABE', 'VRER-ABE', 'ABEsa', 'Sa(KKH)-ABE']

		# name : (pam, 5'or3'pam, gide length, approximated site of DSB site, notes )
		# HDR most effcient within 1-7 bases outside of the DSB, so keeping this will remain standard with non-BE
		self.editor_pamlib = {'spCas9': ('NGG', False, 20, -2, 'requirments work for SpCas9-HF1, eSpCas9 1.1'),
		                      'fnCas9': ('NGG', False, 21, -2, 'highly specifc yet large enzyme'),
		                      'saCas9': ('NNGRRT', False, 21, -2, 'Cas9 S. Aureus 21 base guide'),
		                      'Nme2Cas9c': ('NNNNCC', False, 23, -2, ''),
		                      'spG': ('NGN', False, 20, -2, '20bp-NGN - SpG'),
		                      'SpRY-HighE': ('NRN', False, 20, -2, 'High Efficiency Pam'),
		                      'scCas9': (
			                      'NNGT', False, 20, -2,
			                      '20bp-NNGT - Cas9 S. canis - high efficiency PAM, recommended'),
		                      'stCas9': ('NNAGAA', False, 20, -2, 'Cas9 S. Thermophilus'),
		                      'iSpyMacCas9': ('NAA', False, 20, -2, ''),
		                      'CasX': ('TTCN', True, 20, 18, 'plmCas12e,dltCas12e'),
		                      'AsCas12a': ('TTTV', True, 23, 22, 'TTT(A/C/G)-23bp - Cas12a (Cpf1)'),
		                      'LbCas12a': ('TTTV', True, 23, 22, 'LbCpf1'),
		                      'Cas12d': ('TA', True, 18, 17, 'CasY'),
		                      'Cas12c1': ('TG', True, 23, 22, 'C2c3'),
		                      'AacCas12b': ('TTN', True, 20, 18, ''),
		                      'UnCas12c2': ('TN', True, 20, 18, 'only binds, does not cut')
		                      }

		## ------------Defaults and settings------------------##
		##configure editor options
		self.search_params = self.configure_search_params()

		# configure BE options
		if self.BEmode != 'off':
			self.BE_search_params = self.set_BE_params()
		# ---------------Flags--------------------------#
		self.clininfo_flag = False

		# ---------------Ouputs--------------------------#
		self.all_variant = pd.DataFrame()
		self.all_gene = pd.DataFrame()
		self.all_guides = {}
		self.all_BE = {}

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
			self.search_params = {'spCas9': ('NGG', False, 20, -2, 'SpCas9, SpCas9-HF1, eSpCas9 1.1'),
			                      'saCas9': ('NNGRRT', False, 21, -2, 'Cas9 S. Aureus 21 base guide'),
			                      'CasX': ('TTCN', True, 20, 18, 'plmCas12e,dltCas12e'),
			                      'Cas12a': ('TTTV', True, 23, 22, 'LbCas12,AsCas12a')}

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

	def set_BE_params(self):
		# sets base editor search params, each key is a list of 2 or more; refernce seq search params,
		# then any set that follows starts with the conversion (ex. 'AG' is A --> G) and then the base editors that have the same params

		BE_lib = {'spCas9-def': [('NGG', False, 20, [4, 8]), ('CT', 'BE'), ('AG', 'ABE')],
		          'spCas9-BE14': [('NGG', False, 20, [4, 8]), ('CT', 'BE1|BE3|BE4|HF-BE')],
		          'spCas9-ABE7.9': [('NGG', False, 20, [4, 9]), ('AG', 'ABE7.9')],
		          'spCas9-ABE7.10': [('NGG', False, 20, [4, 7]), ('AG', 'ABE7.10')],
		          'spCas9-max': [('NGG', False, 20, [4, 8]), ('CT', 'BE4max'), ('AG', 'ABEmax')],
		          'Target-AID': [('NGG', False, 20, [2, 8]), ('CT', 'Target-AID')],
		          'spCas9-YE1': [('NGG', False, 20, [4, 7]), ('CT', 'YE1-BE3')],
		          'spCas9-YE': [('NGG', False, 20, [5, 6]), ('CT', 'EE-BE3|YE2-BE3|YEE-BE3')],
		          'spCas9-VQR': [('NGA', False, 20, [4, 8]), ('CT', 'VQR-BE3'), ('AG', 'VQR-ABE')],
		          'spCas9-VRER': [('NGCG', False, 20, [4, 8]), ('CT', 'VRER-BE3'), ('AG', 'VRER-ABE')],
		          'saCas9-BE': [('NNGRRT', False, 21, [3, 12]), ('CT', 'SaBE3', 'SaBE4')],
		          'saCas9-KKh-BE': [('NNNRRT', False, 21, [3, 12]), ('CT', 'Sa(KKH)-BE3')],
		          'saCas9-KKH-ABE': [('NNGRRT', False, 21, [8, 18]), ('AG', 'ABESa', 'Sa(KKH)-ABE')]}

		if self.BEmode == 'default':
			self.BE_search_params = {'spCas9-def': BE_lib['spCas9-def']}

		elif self.BEmode == 'all':
			self.BE_search_params = BE_lib

		else:
			if self.BEmode not in self.BE_choices:
				print('That is not a valid Base Editor')
				print(f'please choose from {self.BE_choices}')
			else:
				for k, v in BE_lib.values():
					if self.BEmode in v[1][-1]:
						self.BE_search_params = {self.BEmode: BE_lib[k][0:2]}

					if len(v) == 3:
						if self.BEmode in v[2][-1]:
							self.BE_search_params = {self.BEmode: BE_lib[k][0] + BE_lib[k][2]}

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
		if 'Doench Score' in df.columns:
			temp = df[df['Editor'] == 'spCas9'].sort_values(by='Doench Score', ascending=False)
			df = pd.concat([temp, df[df['Editor'] != 'spCas9']]).reset_index(drop=True)
		df['Guide_ID'] = [y + str(x) for x, y in zip(list(df.index), list(df['Guide_ID']))]
		# nameout = 'BaseEditors_found.csv' if gtype == 'BE' else 'Guides_found.csv'
		# datenow = date.today().strftime('%Y-%m-%d')
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

		tempgene = pd.read_csv(f"{self.processed_tables}/gene_tables/gene_tables.csv.gz")
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
	def extract_seqs(searchseq, pos, alt, window=30):
		"""
		extracts the sequence +/-30bp surrounding a SNV then swaps ref for alt allele
		"""
		extracted_seq = str(searchseq[pos - window:pos + window])
		extracted_seq = Seq(extracted_seq[0:window] + alt + extracted_seq[window + 1:]).upper()
		return extracted_seq

	def get_refseq_entry(self, term, field):
		'''
		Using ncbiRefSeq.txt to find cds features by either interval, gene name or transcript ID
		example input:
		term, field = 'NM_000532.5', 'tid'
		term, field = 'ENST00000251654.9', 'eid'
		term, field = 'PCCB','name'
		term,field =  'chr3:136250339-136330169','interval'
		'''

		global entry
		labels = ['eid', 'tid', 'chrom', 'strand', 'txStart', 'txEnd',
		          'cdsStart', 'cdsEnd', 'exonStarts', 'exonEnds', 'name', 'exonFrames']

		if field != 'interval':
			not_found = True
			for line in gzip.open(self.annote_path, 'rt'):
				tokens = line.split('\t')
				entry = dict(zip(labels, tokens))
				if term in entry[field]:
					not_found = False
					break

			if not_found:
				entry = None
				print(f"{term} not found in refseq data")

				if '.' in term:
					new_term = term.split('.')[0]
					print(f'searching for {new_term} instead')
					entry = self.get_refseq_entry(new_term, field)

		else:  # only used for intervals search
			not_found = True
			ch = term.split(":")[0]
			start, end = term.split(":")[1].split('-')
			pos = int((int(start) + int(end)) / 2)

			for line in gzip.open(self.annote_path, 'rt'):
				tokens = line.split('\t')
				entry = dict(zip(labels, tokens))
				if ch == entry['chrom']:
					if pos in range(int(entry['txStart']), int(entry['txEnd'])):
						not_found = False
						break
			if not_found:
				entry = None
				print(f"{term} not found in refseq data")

		return entry

	@staticmethod
	def find_codons(dist_from_cds_start, strand):
		'''
		Finds reading frame of SNV in extracted sequence
		'''
		rf = 1 if dist_from_cds_start % 3 == 2 else 2 if dist_from_cds_start % 3 == 0 else 0
		# if strand == '-':
		#	rf = rf * -1
		print(rf, dist_from_cds_start, strand)
		return rf

	@staticmethod
	def get_cds_info(tx_seq, entry):
		'''
		uses entry info to find cds (without utr's)
		'''
		exon_starts = entry['exonStarts'][:-1].split(',')
		exon_ends = entry['exonEnds'][:-1].split(',')
		exon_frames = entry['exonFrames'].replace("\n", "")[:-1].split(',')
		tx_start = int(entry['txStart'])

		exons = [(int(exon_starts[i]) - tx_start, int(exon_ends[i]) - tx_start) for i in range(len(exon_ends))]
		for i in range(len(exon_frames)):
			if exon_frames[i] == '-1':  # -1 means entire exon is UTR
				exons = exons[1:]
				exon_starts = exon_starts[1:]
			else:
				break
		for i in range(1, len(exon_frames)):
			if exon_frames[-i] == '-1':
				exons = exons[0:len(exons) - 1]
				exon_ends = exon_ends[0:-1]
			else:
				break

		# Determine the stop and start of UTR
		exons[0] = (int(entry['cdsStart']) - int(exon_starts[0]) + exons[0][0], exons[0][1])
		exons[-1] = (exons[-1][0], exons[-1][1] - (int(exon_ends[-1]) - int(entry['cdsEnd'])))

		cds = Seq(''.join([str(tx_seq)[a:b] for a, b in exons]))
		if entry['strand'] == '-':
			cds = cds.reverse_complement()

		# translation = cds.translate()
		return [exons, tx_seq, cds]

	def find_transcript_info(self, term, fasta):
		'''
		Using a Refseq Transcript_ID, Ensembl Transcript_ID or coordinates find transcript annotations and transcript sequence
		from either a genome fasta path or given genome sequence
		'''
		# id= 'NM_000532.5' or 'ENST00000251654.9'
		# fasta = f"/groups/clinical/projects/clinical_shared_data/hg38/hg38_chr20.fa.gz"
		if type(fasta) == str:
			fasta_seq = SeqIO.read(gzip.open(fasta, 'rt'), 'fasta')
		else:
			fasta_seq = fasta
		field = 'eid' if term.startswith('E') else 'tid' if term.startswith('N') else 'interval'

		entry = self.get_refseq_entry(term=term, field=field)
		if entry != None:
			tx_seq = fasta_seq.seq[int(entry['txStart']):int(entry['txEnd'])]
			tid_info = self.get_cds_info(tx_seq, entry)
		else:
			tid_info = None
		return entry, tid_info

	def find_snvseq_info(self, snvpos, alt, tid_info, entry, window=30):
		# returns - sequence,feature,translation(if needed)
		# feature: non-coding, utr5,ut3,intron,exon, start_codon, stop_codon
		# snvpos, alt = 11576257, 'T'
		global dist_from_cds_start
		seq, feature, rf = None, None, None

		exons, tx_seq, cds = tid_info
		t_snvpos = int(snvpos) - int(entry['txStart'])
		cdstart, cdsend = int(entry['cdsStart']) - int(entry['txStart']), int(entry['txEnd']) - int(entry['txStart'])
		strand = entry['strand']

		seq = self.extract_seqs(searchseq=tx_seq, pos=t_snvpos - 1, alt=alt, window=30)

		if t_snvpos < 0:
			# not in transcript - shouldn't happen or else no entry would be found
			feature = 'non-coding'

		else:
			if t_snvpos in range(cdstart, cdsend + 1):
				# in CDS
				# find if utr
				feature = 'intron'
				if t_snvpos in range(cdstart, exons[0][0] + 1):
					feature = '3utr' if strand == '-' else '5utr'
				elif t_snvpos in range(exons[-1][1], cdsend + 1):
					feature = '5utr' if strand == '-' else '3utr'

				# find if exon or intron
				else:
					exon_n = 0
					for x in exons:
						# if in exon find reading frame
						if t_snvpos in range(x[0], x[1] + 1):
							# stop and start codon
							feature = 'exon'
							dist = sum([e[1] - e[0] for e in exons[0:exon_n]])
							dist_from_cds_start = dist + (t_snvpos - x[0])

							if strand == '-':
								dist_from_cds_start = (len(cds) - dist_from_cds_start) + 1

							if dist_from_cds_start < 3:
								feature = 'start_codon'

							if dist_from_cds_start > len(cds) - 3:
								feature = 'stop_codon'

							rf = self.find_codons(dist_from_cds_start, strand)
							break
						exon_n += 1
			else:
				# in transcript but not in cds
				seq = str(seq)

				if len(SeqUtils.nt_search(seq[window - 6:window + 5], 'TTTATT')) > 1 or len(
						SeqUtils.nt_search(seq[window - 6:window + 5], 'AATAAA')) > 1:
					feature = 'polya'
				elif len(SeqUtils.nt_search(seq[window - 6:window + 5], 'TATAAA')) > 1 or len(
						SeqUtils.nt_search(seq[window - 6:window + 5], 'ATATTT')) > 1:
					feature = 'promoter'
				elif len(SeqUtils.nt_search(seq[window - 7:], 'GGNCAATCT')) > 1:
					if len(SeqUtils.nt_search(seq[window - 7:window + 6], 'GGNCAATCT')) > 1:
						feature = 'promoter'
					else:
						feature = 'TSS'
				elif len(SeqUtils.nt_search(seq[window + 6:], 'AGATTGNCC')) > 1:
					if len(SeqUtils.nt_search(seq[window - 7: window + 6], 'AGATTGNCC')) > 1:
						feature = 'promoter'
					else:
						feature = 'TSS'
				else:
					feature = 'flanking'

		return seq, feature, rf

	def fetch_query_info(self):
		# Gets Transcript info
		global term
		snv_info = {}

		# If quering by HGVSID with no other info then need to get chromsome/location/alt/ref
		if self.qtype == 'hgvs' and self.hgvscoord is None:
			print("Looking up HGVS in Clinvar.......")

			hgvs_tab = pd.read_csv(self.HGVSlookup_path)
			q_prefixes = [x.split(':')[0] for x in self.queries]
			chroms = set(hgvs_tab.loc[hgvs_tab['TranscriptID'].isin(q_prefixes), 'Chr'])

			for ch in chroms:
				df = pd.read_csv(f"{self.processed_tables}/variant_tables/{ch}_variant.txt")
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
			fasta = SeqIO.read(gzip.open(self.fasta_path.replace('.fa.gz', f'_chr{str(ch)}.fa.gz'), 'rt'), 'fasta')
			new_data = []

			for d in data:
				query, snvpos, ref, alt = d
				print(query, ref, alt, snvpos)
				if self.qtype == 'hgvs':  # pull refseqID from HGVS and search transcript by this
					term = query.split(':')[0]
				if self.qtype == 'coord':  # else use coordsinates to search trancript
					term = f"chr{str(ch)}:{str(snvpos)}-{str(snvpos)}"
				entry, tid_info = self.find_transcript_info(term=term, fasta=fasta)
				if entry != None:
					extracted_seq, feature_annotation, codons = self.find_snvseq_info(snvpos, alt, tid_info, entry,
					                                                                  window=30)
					strand = entry['strand']
					tid, eid = entry['tid'], entry['eid']
				else:
					feature_annotation = 'undetermined/non-coding'
					codons = 'None'
					strand = '+'
					extracted_seq = self.extract_seqs(fasta.seq, snvpos, alt, window=30)
					tid, eid = term, '-'
				new_data.append(
					[query, tid, eid, strand, ref, alt, feature_annotation, extracted_seq, codons,
					 f"chr{str(ch)}:{str(snvpos)}"])

				print('Query term & annoation:', query, feature_annotation)
			snv_info[ch] = new_data

		self.snv_info = snv_info
		# self.write_snv_site_info()

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

	def run_FetchGuides(self, outfile_path):
		global dh, query
		self.fetch_query_info()
		print('Finding Guides.....')
		for ch, data in self.snv_info.items():

			for d in data:
				query, tid, eid, strand, ref, alt, feature_annotation, extracted_seq, codons, coord = d
				dh = DataHandler(query, strand, ref, alt, feature_annotation, extracted_seq, codons, coord)

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

					print(len((guides['gRNA'])), ' guides found for ', query)
			else:
				print(f"No guides found for the query {query}")

		guidedf, BEdf = None, None

		if len(self.all_guides.keys()) != 0:
			guidedf = self.write_guide_csv(self.all_guides, outfile_path)
			self.clininfo_flag = True
			# self.add_clininfo()

		if len(self.all_BE.keys()) != 0:
			BEdf = self.write_guide_csv(self.all_BE, outfile_path)
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
	# Non-dependent
	# === Params ===
	resultsfolder = set_export(str(snakemake.params.main_out))
	gene_report = f"{resultsfolder}/{str(snakemake.params.gene_report)}"
	variant_report = f"{resultsfolder}/{str(snakemake.params.variant_report)}"
	be_report = f"{resultsfolder}/{str(snakemake.params.be_report)}"
	guides_report = f"{resultsfolder}/{str(snakemake.params.guides_report)}"
	# == Intermediate paths
	intermediate_out = set_export(str(snakemake.params.intermediate_out))
	guide_search_params_path = f"{intermediate_out}/{str(snakemake.params.guide_search_params)}"
	snv_site_info_path = f"{intermediate_out}/{str(snakemake.params.snv_site_info)}"
	# == Processed tables branch
	datadir = str(snakemake.params.support_tables)
	annote_path = str(snakemake.params.annote_path)
	# === Wildcards ===
	jobname = str(snakemake.wildcards.job_name)
	# Paths---------------------------
	# input_file = "/groups/clinical/projects/editability/medit_queries/medit_test/test_in/hgvs_test_queries.csv"
	# datadir = "/groups/clinical/projects/editability/tables/"
	# resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"
	# fasta_path = "/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa.gz"
	# annote_path =  "/groups/clinical/projects/editability/tables//processed_tables/ncbiRefSeq.txt.gz"

	# == DEBUG BLOCK ==
	qtype = 'hgvs'
	BEmode = 'off'
	editor = 'all'
	## == == ==
	# == Input Setup ==
	df = pd.read_csv(input_file)
	queries = list(df.iloc[:, 0])
	# == Define guides output path ==
	guides_report_out = be_report if BEmode == 'on' \
		else guides_report

	# == Report processed input variables ==
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
	OUTPUTS TO:
		--> {resultsfolder}
	""")
	# == Get query items ==
	fg = Fetch_Guides(queries,
	                  qtype,
	                  editor,
	                  BEmode,
	                  datadir,
	                  fasta_path,
	                  annote_path
	                  )
	# == Set up object and run core methods ==
	exports = fg.run_FetchGuides(guides_report_out)

	# == Export Intermediate files ==
	fg.write_snv_site_info(snv_site_info_path)
	fg.write_gsearch_params(guide_search_params_path)

	# == Export Variant and Gene tables ==
	fg.add_clininfo(gene_report, variant_report)
# for item_name in exports:
# 	filepath = f"{resultsfolder}/{item_name}.csv"
# 	exports[item_name].to_csv(filepath)


if __name__ == "__main__":
	main()

# Native Modules
import time
import regex as re
import os
from datetime import date
import pandas as pd
import gzip
from Bio.Seq import Seq
from Bio import SeqIO

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
                 editor: str|list,
                 BEmode: str|list,
				 resultsfolder: str,
                 datadir:str,
                 fasta_path:str,
                 annote_path:str,
				 **kwargs):
		"""
			:param queries: list of query terms, either in hgvs format - 'NM_000518.5:c.114G>A' or coords 'chr11:5226778C>T' (COORDS ALLELES MUST BE PLUS STRAND!!)
			:param qtype: 'hgvs' or 'coord'
			  ---> if 'hgvs', providing the coordinates in the kwargs with 'hgvscoord' can reduce processing time
			  ---> hgvs assumes the query is already in clinvar and will generate a variant report with the gene report, if 'coord' then just gene report is created
			:param editor: 'all', 'custom', selected list or str() of the editor choices
			--> custom must contain kwargs, see below
			:param BEmode: 'off','default','all', or select BE editor for base editor choices below
			:param genome: genome used
			:param datadir: folder where tables and pre-computed data live
			:param fasta_path: *Unsure using chromsome seperate files right now but unsure if this will be permenant
			:param kwargs: 'hgvscoord' , 'Jobname','clin_report','gene_report'

			** if 'custom' selcted as editor in kwargs must include pam, pamISFirst, window_size (optional:name)

        """
	##-----------------User Inputs--------------------##
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
		self.job_name = None
		self.kwargs = kwargs

		if 'hgvscoord' in kwargs.keys():
			self.hgvscoord = self.validate_coord(kwargs['hgvscoord'])  #'chr11:5226778C>T'
		if 'gene_report' in kwargs.keys():
			self.gene_report = kwargs['gene_report']
		if 'clin_report' in kwargs.keys():
			self.clin_report= kwargs['clin_report']
		if 'job_name' in kwargs.keys():
			self.job_name = kwargs['job_name']


		#input paths/folders
		self.processed_tables = f"{datadir}/processed_tables/"  # folder with cleaned clinvar/hpa tabs
		self.HGVSlookup_path = f"{self.processed_tables}HGVSlookup.csv"
		self.fasta_path = fasta_path
		self.annote_path = annote_path

		#output paths/folder
		self.resultsfolder = resultsfolder

		#other variables
		self.snv_info = {} # {chrom: (id,snv_pos,ref,alt)}

		##---------------libraries and keys--------------------##
		# [editor]: pam, pamISfirst, win_size, guidelen, scoring, notes/altnames
		self.editor_choices = ['spCas9', 'saCas9', 'spG', 'SpRY-HighE', 'LbCpf1',
							   'scCas9', 'stCas9', 'iSpyMacCas9', 'CasX', 'Cpf1']

		self.BE_choices = ['BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-Gam', 'YE1-BE3', 'EE-BE3', 'YE2-BE3',
						   'YEE-BE3','VQR-BE3','VRER-BE3','SaBE3','SaBE4','SaBE4-Gam','Sa(KKH)-BE3','xBE3','eA3A-BE3',
						   'A3A-BE3','BE-PLUS','ABE7.9','ABE7.10','xABE,ABESa','VQR-ABE','VRER-ABE','Sa(KKH)-ABE']


		# name : (pam, 5'or3'pam, protospace length, approximated site of DSB site, notes )
		# HDR most effcient within 1-7 bases outside of the DSB, so keeping this will remain standard with non-BE
		self.editor_pamlib = {'spCas9': ('NGG', False, 20, -2, 'Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
					  'saCas9': ('NNGRRT', False,21, -2, 'Cas9 S. Aureus 21 base guide'),
					  'spG': ('NGN', False, 20, -2, '20bp-NGN - SpG'),
					  'SpRY-HighE': ('NRN', False,20, -2, 'High Efficiency Pam'),
						'scCas9':('NNGT',False,20,-2,'20bp-NNGT - Cas9 S. canis - high efficiency PAM, recommended'),
					  'stCas9': ('NNAGAA', False,20, -2, 'Cas9 S. Thermophilus'),
					  'iSpyMacCas9': ('NAA', False,20, -2, ''),
					  'CasX': ('TTCN', True, 20, 18, 'Cas12e'),
					  'AsCas12a': ('TTTV', True, 23, 22, 'TTT(A/C/G)-23bp - Cas12a (Cpf1)'),
					  'LbCas12a': ('TTTV', True, 23, 22, 'LbCpf1'),
					  'Cas12c1': ('TG', True, 23, 22, 'C2c3'),
							  }


		## ------------Defaults and settings------------------##
		##configure editor options
		self.search_params = self.configure_search_params()

		#configure BE options
		if self.BEmode != 'off':
			self.BE_search_params = self.set_BE_params()

		#---------------Ouputs--------------------------##
		self.all_variant = pd.DataFrame()
		self.all_gene = pd.DataFrame()
		self.all_guides = {}
		self.all_BE = {}

	def set_guidelen(self, guidelen):
		self.guidelen = guidelen

	def set_win_size(self, win_size):
		self.win_size = win_size

	def configure_search_params(self):
        '''
        set parameteres for the selected editor or editors(not BE editors)
        '''
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


	def set_BE_params(self):
		# sets base editor search params, each key is a list of 2 or more; refernce seq search params,
		# then any set that follows starts with the conversion (ex. 'AG' is A --> G) and then the base editors that have the same params

		BE_lib = {'spCas9-def': [('NGG', False, 20,[4,8]),('CT','BE3', 'BE4', 'BE4max', 'BE4-Gam'),('AG','ABE7.9','ABE7.10','ABEmax')],
				  'spCas9-YE': [('NGG', False, 20, [5, 6]), ('CT', 'YE1-BE3', 'EE-BE3', 'YE2-BE3', 'YEE-BE3')],
				  'spCas9-VQR': [('NGA', False, 20, [4, 8]), ('CT', 'VQR-BE3'), ('AG', 'VQR-ABE')],
				  'spCas9-VRER': [('NGCG', False, 20, [4, 8]), ('CT', 'VRER-BE3'), ('AG', 'VRER-ABE')],
					   'saCas9-KKH': [('NNGRRT', False, 20, [3, 12]),('CT','SaBE3','SaBE4','SaBE4-Gam,Sa(KKH)-BE3'),('AG','ABESa','Sa(KKH)-ABE')]}

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
		if 'Doench Score' in df.columns:
			temp = df[df['Editor'] == 'spCas9'].sort_values(by = 'Doench Score', ascending=False)
			df = pd.concat([temp,df[df['Editor'] != 'spCas9']]).reset_index(drop=True)
		df['Guide_ID'] = [y + str(x) for x,y in zip(list(df.index), list(df['Guide_ID']))]
		nameout = 'BaseEditors_found.csv' if gtype == 'BE' else 'Guides_found.csv'

		datenow = date.today().strftime('%Y-%m-%d')
		if self.job_name != None:
			out = f'{self.resultsfolder}{self.job_name}_{datenow}_{nameout}'
		else:

			out = f'{self.resultsfolder}{datenow}_{nameout}'
		df.to_csv(out,index = False)
		return df

	def add_clininfo(self):
		all_tids = []

		for ch, data in self.snv_info.items():
			all_tids += [d[1] for d in data]

			if self.qtype == 'hgvs':
				tempvar = pd.read_csv(f"{self.processed_tables}/variant_tables/{ch}_variant.txt")
				tempvar = tempvar.loc[tempvar['HGVS_Simple'].isin(list(self.queries))]
				self.all_variant= pd.concat([self.all_variant,tempvar])

		tempgene = pd.read_csv(f"{self.processed_tables}/gene_tables/gene_tables.csv.gz")
		self.all_gene = tempgene.loc[tempgene['TranscriptID'].isin(list(all_tids))]



		datenow = date.today().strftime('%Y-%m-%d')
		self.all_gene.to_csv(f'{self.resultsfolder}/{datenow}_Gene_Report.csv',index = False)

        if self.qtype == 'hgvs':
		self.all_variant.to_csv(f'{self.resultsfolder}/{datenow}_Variant_Report.csv',index = False)



	def extract_seqs(self,searchseq, pos, alt, window=30):
		"""
		extracts the sequence +/-30bp surrounding a SNV then swaps ref for alt allele
		"""
		extracted_seq = str(searchseq[pos - window:pos + window])
		extracted_seq = Seq(extracted_seq[0:window] + alt + extracted_seq[window + 1:]).upper()
		return extracted_seq


	def get_refseq_entry(self,term, field):
		'''
		Using ncbiRefSeq.txt to find cds features by either interval, gene name or transcript ID
		example input:
		term, field = 'NM_000532.5'', 'tid'
		term, field = 'ENST00000251654.9', 'eid'
		term, field = 'PCCB','name'
		term,field =  'chr3:136250339-136330169','interval'
		'''

		labels = ['eid', 'tid', 'chrom', 'strand', 'txStart', 'txEnd',
				  'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
				  'score', 'name', 'cdsStartStat', 'cdsEndStat',
				  'exonFrames']

		if field != 'interval':
			not_found = True
			for line in gzip.open(self.refseq_path, 'rt'):
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

			for line in gzip.open(self.refseq_path, 'rt'):
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


	def find_codons(self,dist_from_cds_start, strand):
		'''
		Finds reading frame of SNV in extracted sequence
		'''
		rf = 1 if dist_from_cds_start % 3 == 2 else 2 if dist_from_cds_start % 3 == 0 else 0
		if strand == '-':
			rf = rf * -1
		return rf


	def get_cds_info(self,tx_seq, entry):
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
		exons[0] = (int(entry['cdsStart']) - int(exon_starts[0]), exons[0][1])
		exons[-1] = (exons[-1][0], exons[-1][1] - (int(exon_ends[-1]) - int(entry['cdsEnd'])))

		cds = Seq(''.join([str(tx_seq)[a:b] for a, b in exons]))
		if entry['strand'] == '-':
			cds = cds.reverse_complement()

		# translation = cds.translate()
		return [exons, tx_seq, cds]


	def find_transcript_info(self,term, fasta):
		'''
		Using a Refseq Transcript_ID, Ensembl Transcript_ID or coordinates find transcript annotations and transcript sequence
		from either a genome fasta path or given genome sequence
		'''
		# id= 'NM_000532.5' or 'ENST00000251654.9'
		# fasta = f"/groups/clinical/projects/clinical_shared_data/hg38/hg38_chr3.fa.gz"
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


	def find_snvseq_info(self,snvpos, alt, tid_info, entry, window=30):
		# returns - sequence,feature,translation(if needed)
		# feature: non-coding, utr5,ut3,intron,exon, start_codon, stop_codon
		# snvpos, alt = 11576257, 'T'
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
							seq_window = (t_snvpos - window, t_snvpos + window)
							dist = sum([e[1] - e[0] for e in exons[0:exon_n]])
							if strand == '+':
								dist_from_cds_start = dist + (t_snvpos - x[0])

							if strand == '-':
								dist = len(cds) - dist
								dist_from_cds_start = dist - (x[1] - t_snvpos)

							if dist_from_cds_start < 3:
								feature = 'start_codon'
							if len(cds) - dist_from_cds_start < 3:
								feature = 'stop_codon'

							rf = self.find_codons(dist_from_cds_start, strand)
							break
						exon_n += 1
			else:
				# in transcript but not in cds
				if seq[window - 6:window + 5].find('TTTATT') > 0 or seq[window - 6:window + 5].find('AATAAA') > 0:
					feature = 'polya'
				elif seq[window - 6:window + 5].find('TATAAA') > 0 or seq[window - 6:window + 5].find('ATATTT') > 0:
					feature = 'promoter'
				elif re.search('GG(A|T|C|G)CAATCT', str(seq[window - 7:])):
					if re.search('GG(A|T|C|G)CAATCT', str(seq[window - 7:window + 6])):
						feature = 'promoter'
					else:
						feature = 'TSS'
				elif re.search('AGATTG(A|T|C|G)CC', str(seq[:window + 6])):
					if re.search('AGATTG(A|T|C|G)CC', str(seq[window - 7:window + 6])):
						feature = 'promoter'
					else:
						feature = 'TSS'
				else:
					feature = 'flanking'

		return seq, feature, rf

	def fetch_query_info(self):
		#Gets Transcript info
		snv_info = {}

		# If quering by HGVSID with no other info then need to get chromsome/location/alt/ref
		if self.qtype == 'hgvs' and self.hgvscoord == None:
			print("Looking up HGVS in Clinvar.......")

			hgvs_tab = pd.read_csv(self.HGVSlookup_path)
			q_prefixes = [x.split(':')[0] for x in self.queries]
			chroms = set(hgvs_tab.loc[hgvs_tab['TranscriptID'].isin(q_prefixes),'Chr'])

			for ch in chroms:
				df = pd.read_csv(f"{self.processed_tables}variant_tables/{ch}_variant.txt")
				gadf = df.loc[df['HGVS_Simple'].isin(self.queries)]
				snv_info[ch] = gadf[['HGVS_Simple', 'PositionVCF', 'RefAlleleVCF', 'AltAlleleVCF']].to_dict('tight')['data']

		#Else All information is given to find transcript info
		else:
			coords = self.queries if self.qtype == 'coord' else self.hgvscoord

			coord_fmt = r'chr[0-9MTXY]*:(\d*)([ATCG]{1})\>([ATCG]{1})'
			for x in range(len(self.queries)):
				ch = coords[x].split(':')[0].replace('chr','')
				if ch not in snv_info.keys():
					snv_info[ch] = []
				snvpos, alt, ref = list(re.search(coord_fmt, coords[x]).groups())
				snv_info[ch].append([self.queries[x], int(snvpos),alt,ref])

		self.snv_info = snv_info
		print("Gathering Variant Genomic Annotation Info.......")

		for ch,data in snv_info.items(): # find transcript info
			fasta = SeqIO.read(gzip.open(self.fasta_path.replace('.fa.gz', f'_chr{str(ch)}.fa.gz'), 'rt'), 'fasta') #<-----How I'm search genome info
			new_data= []

			for d in data:
				query,snvpos, ref, alt = d
				if self.qtype == 'hgvs': #pull refseqID from HGVS and search transcript by this
					term = query.split(':')[0]
				if self.qtype == 'coord': #else use coordsinates to search trancript
					term = f"chr{str(ch)}:{str(snvpos)}-{str(snvpos)}"
				entry, tid_info = self.find_transcript_info(term=term, fasta=fasta)
				extracted_seq, feature_annotation, codons = self.find_snvseq_info(snvpos, alt, tid_info, entry, window=30)
				strand = entry['strand']
				new_data.append([query, entry['tid'],entry['eid'], strand, ref, alt, feature_annotation, extracted_seq, codons, f"chr{str(ch)}:{str(snvpos)}"])
			snv_info[ch] = new_data

		self.snv_info = snv_info


	def validate_hgvs(self,queries):
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

	def validate_coord(self, queries):
		'''
		standardizes input coordinate and checks formatting
		'''
		#q = 'chr11:5226778C>T'
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

	def run_FetchGuides(self):

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
		guidedf, BEdf = None,None
		if len(self.all_guides.keys()) != 0:
		guidedf = self.write_guide_csv(self.all_guides, gtype='guides')
		if len(self.all_BE.keys()) != 0:
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
	# resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"
	# datadir = "/groups/clinical/projects/editability/tables/"
	# fasta_path = "/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa.gz"

	#HGVS TEST
	queries = ['NM_000532.5(PCCB):c.1316A>G (p.Tyr439Cys)', 'NM_000518.5(HBB):c.114G>A', 'NM_000517.6(HBA2):c.99G>A', 'NM_005886.3(KATNB1):c.1A>G']
	qtype = 'hgvs'
	BEmode = 'default'
	editor = 'all'

	# queries = ['chr3:136327650A>G','chr11:5226778C>T','chr16:173128G>A','chr16:57737244A>G']
	# qtype = 'coord'

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

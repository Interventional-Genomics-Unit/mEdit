# Native Modules
import logging
import os
import pickle
import re
# Installed Modules
import pandas as pd
# Project Modules
from dataH import DataHandler
from annotate import Transcript


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
				 editor_request: list,
				 be_request: list,
				 editors: dict,
				 base_editors: dict,
				 dist_from_cutsite: int,
				 datadir: str,
				 fasta_path: str,
				 annote_path: str,
				 **kwargs):
		"""
		:param queries: list of query terms, either in hgvs format - 'NM_000518.5:c.114G>A' or coords 'chr11:5226778C>T'
		(COORDS ALLELES MUST BE PLUS STRAND!!)
		:param qtype: 'hgvs' or 'coord'
		  --> if 'hgvs', providing the coordinates in the kwargs with 'hgvscoord' can reduce processing time
		  --> hgvs assumes the query is already in clinvar and will generate a variant report with the gene report,
		  --> if 'coord' then just gene report is created
		:param editor_request: 'clinical', 'custom', name (list from editor choices)
			--> custom must contain kwargs - pam, pam_is_first,guide_length, dsb_position
		:param be_request: 'off','default',custom, or select BE editor for base editor choices below
			--> custom must contain kwargs - pam, pam_is_first, guide_length, editing_window, target_base,result_base
		:param editors: Dictionary containing information on the current set of editors supported by mEdit
		:param base_editors: Dictionary containing information on the current set of base editors supported by mEdit
		:param dist_from_cutsite: this is how farway the cutsite/dsb_position may be from the start of the variant position
		:param datadir: folder where tables and pre-computed data live
		:param fasta_path: *Unsure using chromsome seperate files right now but unsure if this will be permenant
		:param kwargs:
		"""
		##-----------------User Required_Inputs--------------------##
		if qtype == 'hgvs':
			self.queries = self.validate_hgvs(queries)
		if qtype == 'coord':
			self.queries = self.validate_coord(queries)
		self.qtype = qtype
		self.editor_request = editor_request  ## clinical, custom, list
		self.be_request = be_request  ## off, default, custom, list

		##------------Optional inputs-------------###
		self.dist_from_cutsite, self.window = self.validate_dist_from_cutsite(dist_from_cutsite)

		##----endonuclease and BE search params
		self.editor_lib = editors
		self.be_lib = base_editors
		self.guide_length = 20
		self.pam = str()
		self.pam_is_first = False
		self.dsb_position = int()
		self.editing_window = [4, 8]
		self.conversion = str()
		self.search_params = self.configure_search_params(kwargs)
		self.BE_search_params = self.configure_BE_params(kwargs)

		# input paths/folders
		self.processed_tables = f"{datadir}/processed_tables"  # folder with cleaned clinvar/hpa tabs
		self.HGVSlookup_path = f"{self.processed_tables}/HGVSlookup.csv"
		self.fasta_path = fasta_path
		self.annote_path = annote_path

		##---------------Outputs--------------------##
		self.snv_info = {}  # {chrom: (queryterm,tid,eid,gene,strand,ref,alt,feature,extracted_seq,rf,coord)}
		self.all_guides = {}
		self.all_BE = {}
		self.not_found = {}
		self.all_variant = pd.DataFrame()
		self.all_gene = pd.DataFrame()
		self.nguides_df = pd.DataFrame()

	def configure_search_params(self, kwargs):
		"""
		set paramteres for the selected editor or editors(not BE editors)
		"""
		search_params = {}
		# search for all guides
		if self.editor_request == str('clinical'):
			search_params = self.editor_lib['clinical']
		# set custom editor params
		elif self.editor_request == str('custom'):
			search_params = self.set_custom_params(kwargs)

		# search for selected subset
		else:
			# === At this point, the user-defined list has been thoroughly validated by guide_prediction.py ===
			logging.info(f"EDITOR REQUEST: {self.editor_request} -- Type: {type(self.editor_request)}")
			for editor in list(self.editor_request):
				search_params[editor] = self.editor_lib['all'][editor]

		logging.info(f'Editor(s) set to: {[x for x in search_params.keys()]}')
		return search_params

	def set_custom_params(self, kwargs):
		try:
			self.pam = kwargs['pam'].upper()
		except KeyError:
			logging.error("Custom editor selection MUST contain pam argument")

		try:
			self.pam_is_first = kwargs['pam_is_first']
			if not isinstance(kwargs['pam_is_first'], bool):
				raise TypeError(f"'pam_is_first' must be of type 'bool', got {type(kwargs['pam_is_first'])} instead.")
		except KeyError:
			logging.error("custom editor selection MUST also have  pam_is_first in kwargs")

		try:
			self.guide_length = kwargs['guide_length']
			if not isinstance(kwargs['guide_length'], int):
				raise TypeError(f"'guide_length' must be of type 'int', got {kwargs['guide_length'].__name__} instead.")
		except KeyError:
			logging.error("Custom editor selection MUST also have guide_length in kwargs")
			logging.error(f"Custom editor guide length being set to {self.guide_length}")
			pass
		try:
			dsb_position = str(kwargs['dsb_position']).split(",")
			self.dsb_position = int(dsb_position[0])
		except KeyError:
			logging.error("Custom editor selection MUST also have dsb_loc in kwargs")

		params = {'custom_endonuclease': (self.pam, self.pam_is_first, self.dsb_position, self.guide_length, '')}

		logging.info(
			f"Your custom editor has a {'5 prime PAM' if self.pam_is_first else '3 prime PAM'} set to {self.pam} "
			f"with a spacer length of {self.guide_length} "
			f"5 prime -{self.pam if self.pam_is_first else ''}{'x' * self.guide_length}{'' if self.pam_is_first else self.pam}- 3 prime")
		return params

	def configure_BE_params(self, kwargs):
		# sets base editor search params, each key is a list of 2 or more; refernce seq search params,
		# then any set that follows starts with the conversion (ex. 'AG' is A --> G) and then the base editors that have the same params

		be_search_params = {}
		if self.be_request == 'default':
			be_search_params = self.be_lib['default']
		elif self.be_request == 'off':
			be_search_params = {}
		elif self.be_request == 'custom':
			be_search_params = self.set_be_custom_params(kwargs)
		else:
			# === At this point, the user-defined list has been thoroughly validated by guide_prediction.py ===
			for base_editor in list(self.be_request):
				be_search_params[base_editor] = self.be_lib['all'][base_editor]

		logging.info(f'Base Editor(s) set to: {[x for x in be_search_params.keys()]}')
		return be_search_params

	def set_be_custom_params(self, kwargs):
		try:
			self.pam = kwargs['pam'].upper()
		except KeyError:
			logging.error("Custom editor selection MUST contain pam argument")

		try:
			self.pam_is_first = kwargs['pam_is_first']
			if not isinstance(kwargs['pam_is_first'], bool):
				raise TypeError(f"'pam_is_first' must be of type 'bool', got {type(kwargs['pam_is_first'])} instead.")
		except KeyError:
			logging.error("Custom editor selection MUST also have  pam_is_first in kwargs")
			logging.error(f"Custom editor pam_is_first is being set to {self.pam_is_first}")
			pass

		try:
			self.guide_length = int(kwargs['guide_length'])
		except KeyError:
			logging.error("Custom editor selection MUST also have guide_length in kwargs")
			logging.error(f"Custom editor pam_is_first is being set to {self.guide_length}")
			pass
		try:
			target_base = kwargs['target_base'].upper()
			result_base = kwargs['result_base'].upper()
			if target_base not in ['A', 'G', 'C', 'T'] or target_base not in ['A', 'G', 'C', 'T']:
				raise TypeError("target_base or result_base need to equal 'A','C','G' or 'T'")
			self.conversion = target_base + result_base
		except KeyError:
			logging.error("Custom editor selection MUST also have 'target_base' and 'result_base' in kwargs")

		params = {'custom_be': [(self.pam, self.pam_is_first, self.guide_length, self.editing_window, ''),
								(self.conversion, 'custom_be')]}

		logging.info(
			f"Your custom base editor has a {'5 prime PAM' if self.pam_is_first else '3 prime PAM'} set to {self.pam} "
			f"with a spacer length of {self.guide_length} "
			f"The base editor window is between position {self.editing_window[0]} and {self.editing_window[1]} "
			f"and will convert {self.conversion[0]}---->{self.conversion[1]} "
			f"5 prime -{self.pam if self.pam_is_first else ''}{'x' * self.guide_length}{'' if self.pam_is_first else self.pam}- 3prime")
		return params

	def validate_dist_from_cutsite(self, dist_from_cutsite):
		'''
		Makes sure the extracted genome window is big enough
		to accomodate distance from cutsite
		'''
		if dist_from_cutsite <= 7:
			return dist_from_cutsite, 50
		elif dist_from_cutsite > 7 and dist_from_cutsite <= 200:
			return dist_from_cutsite, 50 + (dist_from_cutsite - 7)
		else:
			logging.warning('The dist_from_cutsite must be under 200'
							'The current dist_from_cutsite = {self.dist_from_cutsite} is being changed  to 200')
			return 200, 50 + (200 - 7)

	def write_besearch_params(self, outfile):
		# writes pickle of selected be guide search params for later use in process_genome
		with open(outfile, 'ab') as gfile:
			pickle.dump(self.BE_search_params, gfile)

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

	def write_not_found(self, outfile):
		pd.Series(data=self.not_found.values(), index=self.not_found.keys()).to_csv(outfile)

	def write_guide_csv(self, guides, outfile):
		df = pd.DataFrame(guides)
		df['Guide_ID'] = [y + str(x) for x, y in zip(list(df.index), list(df['Guide_ID']))]
		df.to_csv(outfile, index=False)
		return df

	def add_clinvar(self, gadf):
		self.all_variant = pd.concat([self.all_variant, gadf])

	def write_reports(self, gene_out, variant_out, nguides_out):
		try:
			guides_found = True if len(self.all_guides['QueryTerm']) > 0 else False
		except KeyError:
			guides_found = False
		# writes
		all_tids = []
		searched_queries = []
		for chrom in self.snv_info.keys():
			for v in self.snv_info[chrom]:
				all_tids.append(v[1])
				searched_queries.append(v[0])

		# creates gene table from queries
		tempgene = pd.read_csv(f"{self.processed_tables}/gene_tables/gene_tables.csv.gz")
		self.all_gene = tempgene.loc[tempgene['TranscriptID'].isin(list(all_tids))]
		self.all_gene.to_csv(gene_out, index=False)

		# tallies number of guides per query
		if guides_found:
			allgdf = pd.DataFrame(self.all_guides)
			nguides_per_query = allgdf.groupby(['QueryTerm']).size().reset_index()
			nguides_per_query = nguides_per_query.rename(columns={0: 'Number of Guides Found'})
			self.nguides_df = pd.DataFrame(searched_queries, columns=['QueryTerm']).join(
				nguides_per_query.set_index('QueryTerm'), on='QueryTerm')
			self.nguides_df['Number of Guides Found'] = self.nguides_df['Number of Guides Found'].fillna(0).astype(
				'int')
			self.nguides_df.to_csv(nguides_out, index=False)

		if self.qtype == 'hgvs':
			self.all_variant.to_csv(variant_out, index=False)

	def add_not_found(self, nfqueries, reason):
		# adds queries that are not found ex: no guides or hgvs not found
		for nf in nfqueries:
			self.not_found[nf] = reason

	@staticmethod
	def extract_seqs(searchseq, pos, ref, alt, window):
		# extracts the sequence +/windowbp surrounding a SNV then swaps ref for alt allele

		alt = alt.lower()
		searchseq = searchseq.upper()

		if len(ref) == len(alt):  ##substitution
			extracted_seq = str(searchseq[pos - window:pos + window])
			variant_seq = (extracted_seq[0:window] + alt + extracted_seq[window + len(alt):])  # .upper()


		elif len(ref) > len(alt):  # deletion
			diff = len(ref) - len(alt)
			extracted_seq = str(searchseq[pos - window:(pos + window + diff)])
			variant_seq = (extracted_seq[0:window] + alt + extracted_seq[window + len(ref):])  # .upper()

		elif len(ref) < len(alt):  # insertion
			extracted_seq = str(searchseq[pos - window:pos + window - (len(alt))])
			variant_seq = (extracted_seq[0:window] + alt + extracted_seq[window:])  # .upper()
		# print(new_seq)

		else:
			logging.warning('Variant in reference and alternative genomes do not comply')
			logging.warning(ref, alt, pos)
		return variant_seq

	def get_chroms(self, queries):
		# finds unique chromosome for each hgvs
		# needs to happen in order to select right fasta file
		hgvs_tab = pd.read_csv(self.HGVSlookup_path)
		q_prefixes = [x.split(':')[0].split('.')[0] for x in queries]
		chroms = set(hgvs_tab.loc[hgvs_tab['HGVS_ID'].isin(q_prefixes), 'Chromosome'])
		return chroms

	def fetch_query_info(self):
		# Gets Transcript info
		# If quering by HGVSID with no other info then need to get chromsome/location/alt/ref
		if self.qtype == 'hgvs':
			logging.info("Looking up HGVS in Clinvar")
			chroms = self.get_chroms(self.queries)

			for ch in chroms:
				df = pd.read_csv(f"{self.processed_tables}/variant_tables/{ch}_variant.txt",
								 low_memory=False)
				gadf = df.loc[df['HGVS_Simple'].isin(self.queries)]
				self.add_clinvar(gadf)
				self.snv_info[ch] = \
				gadf[['HGVS_Simple', 'PositionVCF', 'RefAlleleVCF', 'AltAlleleVCF']].to_dict('tight')['data']

			# record not found
			found = []
			for k in self.snv_info.keys():
				for v in self.snv_info[k]:
					found.append(v[0])
			notfound = list(set(self.queries).difference(set(found)))
			self.add_not_found(notfound, 'hgvs not found in medit variant database')

		# Else All information is given to find transcript info
		else:
			coord_fmt = r'chr[0-9MTXY]*:(\d*)([ATCG]{1,10})\>([ATCG]{1,10})'
			for query in self.queries:
				ch = query.split(':')[0].replace('chr', '')
				if ch not in self.snv_info.keys():
					self.snv_info[ch] = []
				snvpos, alt, ref = list(re.search(coord_fmt, query).groups())
				self.snv_info[ch].append([query, int(snvpos), alt, ref])

	def find_transcript_info(self):

		logging.info("Gathering variant genomic info from RefSeq file")

		for ch, data in self.snv_info.items():  # find transcript info
			new_data = []
			snvcoords = [f'chr{ch}:{d[1]}-{d[1]}' for d in data]
			try:
				Transcript.load_transcripts(self.annote_path, snvcoords)
			except IndexError:
				# TODO: Something on commit f6f4b19 has changed the communication between fetchGuides and annotate.py
				# TODO: @ACTG802 --> please advise on whether this needs fixing or just a bypass message would be enough
				logging.warning(f"WARNING: function 'load_transcripts' from annotate.py got the wrong number of indices from snvcoords"
								f"This is how it looks like: {snvcoords}")
				continue

			# === Transitioning SeqIO.read to direct import of Pickled SeqRecord Objects ===
			chr_fasta_path = f"{self.fasta_path}/chr{str(ch)}.pkl"
			try:
				logging.info(f"Finding transcripts information on chromosome chr{str(ch)}")
				with open(chr_fasta_path, 'rb') as pfile:
					fasta = pickle.load(pfile)
			except FileNotFoundError:
				logging.warning(
					f"The file {chr_fasta_path} was not found. Please regenerate background data and check the target directory")
				continue
			except pickle.UnpicklingError:
				logging.warning(f"The file {chr_fasta_path} is not in the correct format. Please regenerate background data")
				continue

			for d in data:

				query, snvpos, ref, alt = d
				if re.search('[^ATCG]', ref + alt) is None:

					snvcoord = f'chr{ch}:{d[1]}-{d[1]}'
					tx = Transcript.transcript(snvcoord)

					if tx == 'intergenic':
						tid, eid, gname = '-', '-', '-'
						feature_annotation = tx
						rf, strand = 'None', '+'
						extracted_seq = self.extract_seqs(searchseq=fasta.seq, pos=snvpos, ref=ref, alt=alt,
														  window=self.window)


					else:
						eid, tid, gname, strand, txstart = tx.tx_info()
						tx_seq = tx.get_tx_seq(fasta)
						t_snvpos = int(snvpos) - int(txstart)
						extracted_seq = self.extract_seqs(searchseq=tx_seq, pos=t_snvpos - 1, ref=ref, alt=alt,
														  window=self.window)

						if len(extracted_seq) != self.window * 2:
							# if flanking or utr the extracted seq needs to come from the chromosome file
							logging.info('Extracted Sequence Length OUT OF BOUNDS in Transcript Sequence ')
							logging.info('Searching chromosome sequence instead')

							extracted_seq = self.extract_seqs(
								searchseq=fasta.seq[snvpos - self.window * 3:snvpos + self.window * 3],
								pos=self.window * 3,
								ref=ref,
								alt=alt,
								window=self.window)

						feature_annotation, rf = tx.feature, tx.rf

					new_data.append(
						[query, tid, eid, gname, strand, ref, alt, feature_annotation, extracted_seq, rf,
						 f"chr{str(ch)}:{str(snvpos)}"])

				else:
					self.add_not_found(query, 'ref or alt allele contain non-ATCG characters')
			self.snv_info[ch] = new_data

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
		n = len(validated_queries)
		logging.info(f'{n} out of {len(queries)} HGVS IDs validated')
		if n == 0:
			logging.warning('Query are not in the correct HGVS Format')
		return validated_queries

	@staticmethod
	def validate_coord(queries):
		'''
		standardizes input coordinate and checks formatting
		'''
		# q = 'chr11:5226778C>T'
		# coord_fmt = r'(chr[0-9XYM]*:\d*(A|T|C|G)>(A|T|C|G))'
		coord_fmt = r'(chr[0-9XYM]*:\d*([ATCG]{1,10})>([ATCG]{1,10}))'
		validated_queries = []
		for q in set(queries):
			if re.search(coord_fmt, q):
				validated_queries.append(re.search(coord_fmt, q).groups()[0])
			else:
				pass
		n = len(validated_queries)
		logging.info(f'{n} out of {len(queries)} Coordinate IDs validated')
		if n == 0:
			logging.warning('Queries are not in the correct format: Coordinate + allele')
		return validated_queries

	def run_FetchGuides(self, outfile_guides, outfile_be_guides, models_dir):
		# global dh, query
		self.fetch_query_info()
		self.find_transcript_info()
		logging.info('Guide scanning routine initiated')
		for ch, data in self.snv_info.items():

			for d in data:
				try:
					query, tid, eid, gname, strand, ref, alt, feature_annotation, extracted_seq, codons, coord = d
				except ValueError:
					logging.warning(
						f"WARNING: fetchGuides.py --> The query below has the wrong number of values to unpack. Needs further investigation:\n{d}")
					continue
				dh = DataHandler(query, strand, ref, alt, feature_annotation, models_dir, extracted_seq, codons, coord,
								 gname, self.dist_from_cutsite)

				guides, beguides = dh.get_Guides(self.search_params, self.BE_search_params)

				logging.info(f"{len(guides['gRNA']) + len(beguides['gRNA'])} found for {query}")

				if len(guides['gRNA']) > 0:
					if len(self.all_guides.keys()) == 0:
						self.all_guides = guides
					else:
						for k, v in guides.items():
							self.all_guides[k] += v

				if len(beguides['gRNA']) > 0:
					if len(self.all_BE.keys()) == 0:
						self.all_BE = beguides
					else:
						for k, v in beguides.items():
							self.all_BE[k] += v

				if len(guides['gRNA']) + len(beguides['gRNA']) == 0:
					self.add_not_found([query], 'no guides found')

		if len(self.all_guides.keys()) != 0:
			self.write_guide_csv(self.all_guides, outfile_guides)
		else:
			logging.info('No Endonuclease Guides found for any queries')
			print('No Endonuclease Guides found for any queries')
			exit(0)

		if len(self.all_BE.keys()) != 0:
			self.write_guide_csv(self.all_BE, outfile_be_guides)
		else:
			if self.be_request != 'off':
				logging.info('No Base Editor Guides found for any queries')


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	input_file = str(snakemake.input.query_manifest)
	fasta_path = str(snakemake.input.assembly_dir_path)
	# === Outputs ===
	guides_report = str(snakemake.output.guides_report_out)
	gene_report = str(snakemake.output.gene_report)
	variant_report = str(snakemake.output.variant_report)
	nguides_report = str(snakemake.output.nguides_report)
	be_report = str(snakemake.output.be_report_out)
	guide_search_params_path = str(snakemake.output.guide_search_params)
	guide_be_search_params_path = str(snakemake.output.guide_be_search_params)
	snv_site_info_path = str(snakemake.output.snv_site_info)
	guides_not_found_path = str(snakemake.output.guides_not_found_out)
	#   == Log File ==
	logfile_path = str(snakemake.output.logfile_path)
	# === Params ===
	resultsfolder = set_export(str(snakemake.params.main_out))
	#   == Processed tables branch ==
	datadir = str(snakemake.params.support_tables)
	annote_path = str(snakemake.params.annote_path)
	#   == Editor Parameters
	editors_path = str(snakemake.params.editors)
	base_editors_path = str(snakemake.params.base_editors)
	models_path = str(snakemake.params.models_path)
	distance_from_cutsite = int(snakemake.params.distance_from_cutsite)
	#   == Run Parameters ==
	qtype = str(snakemake.params.qtype)
	be_request = snakemake.params.be_request
	editor_request = snakemake.params.editor_request
	#   == Custom Editor Parameters ==
	pam = str(snakemake.params.pam)
	guide_length = int(snakemake.params.guide_length)
	pam_is_first = bool(snakemake.params.pam_is_first)
	dsb_position = int(snakemake.params.dsb_position)
	editing_window = tuple(snakemake.params.editing_window)
	target_base = str(snakemake.params.target_base)
	result_base = str(snakemake.params.result_base)

	# === Log Process Initialization
	# Configure the logging system
	logging.basicConfig(
		level=logging.DEBUG,  # Set the minimum log level (DEBUG logs everything)
		format='%(asctime)s [%(levelname)s] %(message)s',  # Define log format
		handlers=[
			logging.FileHandler(logfile_path),  # Log to file
		]
	)

	logging.info('=== INITIALIZING GUIDE PREDICTION ROUTINE ===')

	# == Create dummy files for optional outputs
	optional_outputs = [be_report, variant_report, gene_report, nguides_report, guides_report]
	for report in optional_outputs:
		# Create an empty DataFrame
		df = pd.DataFrame()
		# Export the empty DataFrame to a CSV file
		df.to_csv(str(report))

	# == Input Setup ==
	df = pd.read_csv(input_file)
	queries = list(df.iloc[:, 0])

	# == Editors / BEs Setup ==
	with open(editors_path, 'rb') as edfile:
		editors = pickle.load(edfile)
	with open(base_editors_path, 'rb') as befile:
		base_editors = pickle.load(befile)

	# == Report processed input variables ==
	# print(f"""
	# Currently running fetchGuides.py
	# INPUT VARIABLES:
	# 	n of Queries: {len(queries)}
	# 	Query Type: {qtype}
	# 	be_request: {be_request}
	# 	editor_request: {editor_request}
	# PATH TO REFERENCE:
	# 	-> {fasta_path}
	# SUPPORT DATA DIRECTORY:
	# 	-> {datadir}
	# OUTPUTS TO:
	# 	--> {resultsfolder}
	# """)

	# == Get query items ==
	fg = Fetch_Guides(queries,
					  qtype,
					  editor_request,
					  be_request,
					  editors,
					  base_editors,
					  distance_from_cutsite,
					  datadir,
					  fasta_path,
					  annote_path,
					  pam=pam,
					  pam_is_first=pam_is_first,
					  guide_length=guide_length,
					  dsb_position=dsb_position,
					  editing_window=editing_window,
					  target_base=target_base,
					  result_base=result_base
					  )
	# == Set up object and run core methods ==
	fg.run_FetchGuides(guides_report, be_report, models_path)

	# == Export Intermediate files ==
	fg.write_snv_site_info(snv_site_info_path)
	fg.write_gsearch_params(guide_search_params_path)
	fg.write_besearch_params(guide_be_search_params_path)

	# == Export Variant, Guide Totals and Gene tables ==
	fg.write_reports(gene_report, variant_report, nguides_report)

	# == Export Not Found table
	fg.write_not_found(guides_not_found_path)

	logging.info('=== GUIDE PREDICTION ROUTINE FINALIZED ===')

if __name__ == "__main__":
	main()

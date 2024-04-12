# Native modules
import re
from typing import Literal
import time
from argparse import ArgumentParser as argp
# Installed modules
import pandas as pd
from Bio import Entrez
import urllib.request
import urllib.error
import urllib3
import xmltodict
# Literals
Orient = Literal["columns", "index", "tight"]


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='ncbi_cross_database',
		description='',
		usage='%(prog)s [options]')
	# Define arguments
	parser.add_argument('input_list',
						help='Path to a plain txt containing a column of NCBI ID entries')  # positional argument
	parser.add_argument('login_email',
						help='Email used to connect with NCBI Entrez')  # positional argument
	parser.add_argument('-o', default='hgvs_list.txt',
						dest='outfile',
						help='Path to output the list of HGVSs')
	parser.add_argument('-f',
						dest='db_from',
						default='gene',
						choices=['clinvar', 'pubmed', 'gene', 'snp'],
						help='The NCBI database where the IDs in the input list will be matched against')
	parser.add_argument('-t',
						dest='db_to',
						default='clinvar',
						choices=['clinvar', 'pubmed', 'gene', 'snp'],
						help='The target NCBI database where the IDs in the input list will be linked to')

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def elink_routine(dbfrom, dbto, hit_uid):
	dup_check = []
	not_found = ""
	linked = ""
	link_record = ""
	handle = None
	server_attempts = 0
	try:
		handle = Entrez.elink(dbfrom=f"{dbfrom}", db=dbto, id=f"{hit_uid}", idtype="uid")
	except urllib.error.HTTPError as err:
		if err.code == 500:
			print(f'An internal server error occurred while handling the accession {hit_uid}')
			not_found = hit_uid
			return linked, hit_uid, not_found
	try:
		link_record = Entrez.read(handle)
	except RuntimeError:
		not_found = hit_uid
	if link_record:
		try:
			linked = []
			# linked = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
			for link_idx in range(len(link_record[0]['LinkSetDb'][0]['Link'])):
				link_id = link_record[0]['LinkSetDb'][0]['Link'][link_idx]['Id']
				linked.append(link_id)
				if link_id not in dup_check:
					dup_check.append(link_id)
		except (IndexError, KeyError):
			not_found = hit_uid
	handle.close()
	return linked, hit_uid, not_found, link_record


def cross_db_search(query_list, db_from, dbto):
	progress = 0
	source2target = {}
	not_found_list = []
	full_record = ''
	for hit in query_list:
		progress += 1
		dup_check = []
		uid = ''
		max_retries = 50
		search_record = {}
		# Standardize protein identifiers to NCBI UIDs through ESearch
		# Introduce a delay of 1 second before making the request

		# Attempt the search with retries
		for _ in range(max_retries):
			try:
				handle = Entrez.esearch(db=f"{db_from}", term=f"{hit}", idtype="acc")
				search_record = Entrez.read(handle)
				print("Processed Esearch request")
				try:
					uid = search_record['IdList'][0]
					print(f"Entry {hit}: Found UID {uid}")
				except IndexError:
					continue
				handle.close()
				break  # Break the loop if successful
			except urllib.error.HTTPError as e:
				if e.code == 429:  # HTTP 429: Too Many Requests
					print(f"Received HTTP 429 error. Retrying in 10 seconds...")
					time.sleep(10)
				else:
					continue  # Re-raise other HTTP errors
			except RuntimeError:
				continue

		# # Standardize identifiers to NCBI UIDs through ESearch
		# handle = Entrez.esearch(db=f"{db_from}", term=f"{hit}", idtype="acc")
		# search_record = Entrez.read(handle)
		# try:
		# 	uid = search_record['IdList'][0]
		# except IndexError:
		# 	continue
		# handle.close()

		# Loop through databases (found in config) and grab Nuccore UIDs
		if uid in set(dup_check):
			continue
		print("Process Elink routine")
		link_list, loop_nuc_acc, not_found_hit, full_record = elink_routine(db_from, dbto, uid)
		if not_found_hit:
			not_found_list.append(not_found_hit)
			continue
		if link_list:
			dup_check.append(uid)
			source2target.setdefault(loop_nuc_acc, link_list)

	return source2target, list(set(not_found_list)), full_record


def url2xml_dict(url):
	http = urllib3.PoolManager()
	file = http.request('GET', str(url))
	data = file.data
	data = xmltodict.parse(data)
	return data


def fetch_hgvs_list(query_list, db_name):
	# Get Genbank records for each Nuccore UID
	hgvs_records = {}
	not_found = []
	for uid in query_list:
		print(f"Processing {uid}")
		handle = Entrez.esummary(db=db_name, id=f"{uid}", rettype="clinvarset", retmode="default")
		xml_data = url2xml_dict(handle.url)
		description = ''
		germline_description = ''
		try:
			if xml_data['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['clinical_impact_classification'][
				'description']:
				description = \
				xml_data['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['clinical_impact_classification'][
					'description']
			if xml_data['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['germline_classification'][
				'description']:
				germline_description = \
				xml_data['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['germline_classification'][
					'description']
			if re.search('pathogenic', description, re.IGNORECASE) or re.search('pathogenic', germline_description,
																				re.IGNORECASE):
				hgvs = xml_data['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['title']
				# Returns a list  of Genbank SeqRecords objects
				hgvs_records.setdefault(uid, hgvs)
				print(f"HGVS found {hgvs}; Description: {description}")
		except KeyError as e:
			print(f"No description or no pathogenic flag found for {uid}")
			not_found.append(uid)
			continue
	# Returns a list  of Genbank SeqRecords objects
	return hgvs_records, not_found


def create_dataframe(data: dict, orient: Orient) -> pd.DataFrame:
	return pd.DataFrame.from_dict(data, orient=orient)


def cross_db():
	args = parse_arguments()
	in_path = args.input_list
	out_path = args.outfile
	db_from = args.db_from
	dbto = args.db_to
	entrez_login = args.login_email
	queries_df = pd.read_csv(in_path, low_memory=False, header=None)
	queries = queries_df.iloc[:, 0].tolist()
	final_df_orient = Orient.__args__[0]

	# === Entrez authentication ===
	Entrez.email = entrez_login

	# === Query NCBI to get UIDs associated with the input IDs using ESearch/ELink ===
	print("Linking query IDs to Nuccore entries")
	hit_to_link, hits_not_found, record = cross_db_search(queries, db_from, dbto)

	# === Translate UID to HGVS if the program is running on default params ===
	print("Cross database pre-processing finalized")
	if dbto == 'clinvar':
		print("Convert input data to HGVS")
		hgvs_list = []
		for input_id in hit_to_link:
			hgvs_list.extend(fetch_hgvs_list(hit_to_link[input_id], dbto))
		hit_to_link = hgvs_list[0]
		final_df_orient = Orient.__args__[1]

	# === Format output table ===
	print("Format output table")
	df_linked = create_dataframe(hit_to_link, orient=final_df_orient)
	df_linked.to_csv(out_path, index=False, header=False)


if __name__ == "__main__":
	cross_db()

# Native modules
import re
from argparse import ArgumentParser as argp
# Installed modules
import pandas as pd
from Bio import Entrez
import yaml
import urllib.request
import urllib.error
import urllib3
import xmltodict


def parse_arguments():
	#  Launch argparse parser
	parser = argp(
		prog='omim2hgvs',
		description='',
		usage='%(prog)s [options]')
	# Define arguments
	parser.add_argument('omim_list',
	                    help='Path to a plain txt containing a column of omim entries')  # positional argument

	# Parse arguments from the command line
	arguments = parser.parse_args()
	return arguments


def elink_routine(db, hit_uid):
	dup_check = []
	not_found = ""
	linked = ""
	link_record = ""
	server_attempts = 0
	try:
		handle = Entrez.elink(dbfrom="omim", db=db, id=f"{hit_uid}")
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


def cross_db_search(query_list, db_list):
	progress = 0
	source2target = {}
	not_found_list = []
	full_record = ''
	for hit in query_list:
		progress += 1
		dup_check = []
		# Standardize protein identifiers to NCBI UIDs through ESearch
		handle = Entrez.esearch(db="omim", term=f"{hit}", idtype="acc")
		search_record = Entrez.read(handle)
		try:
			uid = search_record['IdList'][0]
		except IndexError:
			continue
		handle.close()

		# Loop through databases (found in config) and grab Nuccore UIDs
		for db_name in db_list:
			if uid in set(dup_check):
				continue
			link_list, loop_nuc_acc, not_found_hit, full_record = elink_routine(db_name, uid)
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


def nuc_to_gb(query_list, db_name):
	# Get Genbank records for each Nuccore UID
	gb_records = {}
	for uid in query_list:
		handle = Entrez.esummary(db=db_name, id=f"{uid}", rettype="clinvarset", retmode="default")
		xml_data = url2xml_dict(handle.url)
		try:
			description = xml_data['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['clinical_significance']['description']
			if re.search('pathogenic', description, re.IGNORECASE):
				hgvs = xml_data['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['title']
				# Returns a list  of Genbank SeqRecords objects
				gb_records.setdefault(uid, hgvs)
		except KeyError:
			gb_records.setdefault(uid, xml_data)
			continue
		# Returns a list  of Genbank SeqRecords objects
	return gb_records


def main():

	args = parse_arguments()
	in_path = args.omim_list
	queries_df = pd.read_csv(in_path, low_memory=False, header=None)
	queries = queries_df.iloc[:, 0].tolist()
	# DEBUG
	# queries = ["264300", "610006", "616034", "246450", "605911"]

	# Load config file
	with open("/groups/doudna/projects/daniel_projects/editability/config/id_convert.yaml", "r") as f:
		config = yaml.load(f, Loader=yaml.FullLoader)

	# NCBI databases to search
	efecth_db_list = config['efetch_db']
	# Entrez authentication
	print("Entrez login")
	Entrez.email = config["entrez_login"]

	# Query NCBI to get nuccore UIDs associated with the protein hits using ESearch/ELink
	print("Linking protein hit ids to Nuccore entries")
	hit_to_link, hits_not_found, record = cross_db_search(queries, efecth_db_list)

	# for key in hit_to_link:
	a = nuc_to_gb(hit_to_link['264300'], 'clinvar')


if __name__ == "__main__":
	main()

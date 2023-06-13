# Native modules
import os
import copy
import re
# Installed modules
import urllib.request
import urllib.error
import urllib3
import xmltodict
import yaml
# Biopython
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Load config file
with open("config/entrez.yaml", "r") as f:
	config = yaml.load(f, Loader=yaml.FullLoader)


def get_clinvar_uid(word_query):
	handle = Entrez.esearch(db="clinvar", term=f"{word_query}", idtype="acc")
	search_record = Entrez.read(handle)
	try:
		uid = search_record['IdList'][0]
	except IndexError:
		raise "No HGVSs were found for this search"
	handle.close()
	return uid


def elink_routine(dbfrom, dbto, hit_uid):
	dup_check = []
	not_found = ""
	linked = ""
	link_record = ""
	server_attempts = 0
	try:
		handle = Entrez.elink(dbfrom=dbfrom, db=dbto, id=f"{hit_uid}")
	except urllib.error.HTTPError as err:
		if err.code == 500:
			print(f"An internal server error occurred while handling the accession {hit_uid}")
			not_found = hit_uid
			return linked, hit_uid, not_found
	try:
		link_record = Entrez.read(handle)
	except RuntimeError:
		not_found = hit_uid
	if link_record:
		try:
			linked = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
			if linked not in dup_check:
				dup_check.append(linked)
		except (IndexError, KeyError):
			not_found = hit_uid
	handle.close()
	return linked, hit_uid, not_found


def url2xml_dict(url):
	file = urllib3.request("GET", url)
	data = file.data
	data = xmltodict.parse(data)
	return data


def get_gene_xml(query):
	# Get XML records from Entrez' Gene database
	handle = Entrez.efetch(db="gene", id=f"{query}", rettype="xml", retmode="text")
	xml_data = url2xml_dict(handle.url)
		# Returns a list  of Genbank SeqRecords objects
	return xml_data


def get_gene_coords(gene_xml_dict):
	gene_coords = {}
	chromosomes = gene_xml_dict['Entrezgene-Set']['Entrezgene']['Entrezgene_locus']['Gene-commentary']
	chromosomes_list_length = len(gene_xml_dict['Entrezgene-Set']['Entrezgene']['Entrezgene_locus']['Gene-commentary'])

	# Get gene coordinates for each individual chromosome
	for chr_index in range(0, chromosomes_list_length):
		chromosome_acc = chromosomes[chr_index]['Gene-commentary_accession']
		chromosome_version = chromosomes[chr_index]['Gene-commentary_version']
		chromosome_id = f"{chromosome_acc}.{chromosome_version}"
		gene_coords.setdefault(chromosome_id, [
			gene_xml_dict['Entrezgene-Set']['Entrezgene']['Entrezgene_locus']['Gene-commentary'][chr_index]['Gene-commentary_seqs']['Seq-loc']['Seq-loc_int']['Seq-interval']['Seq-interval_from'],
			gene_xml_dict['Entrezgene-Set']['Entrezgene']['Entrezgene_locus']['Gene-commentary'][chr_index]['Gene-commentary_seqs']['Seq-loc']['Seq-loc_int']['Seq-interval']['Seq-interval_to']
		])
	return gene_coords

# TODO: Editar funcao pra aglutinar todos os GENE UIDs ligados aos seus respectivos SeqRecords -> facilitar busca na funcao get_genbank
# def parse_features(seqrecords_list):


def get_genbank(gene_coords, win_size):
	prot_dict = {}
	record_target = {}
	for chromosome_id in gene_coords:
		# Fetch Genbank entry
		handle = Entrez.efetch(db="nucleotide", id=str(chromosome_id), rettype="gbwithparts", retmode="text")
		record = SeqIO.read(handle, "genbank")
		# Select the 'features' object from the GB file
		features = record.features[0]


		#
		# prot_id = qualifiers["protein_id"][0]
		# # Search for the protein-ids of interest
		# if re.search(prot_id, uid_to_acc[query][hit_uid][0]):

		# Process feature information for future ref
		f_start = int(gene_coords[chromosome_id][0])
		f_end = int(gene_coords[chromosome_id][1])
		f_strand = features.strand
		highlight_feature = copy.deepcopy(features)
		highlight_feature.type = "highlight"
		# Set start/end coords using window size
		start = max(int(min([f_start, f_end])) - win_size, 0)
		end = min(int(max([f_start, f_end])) + win_size + 1, len(record.seq))
		f_len = end - start
	
		# Create a SeqRecord object with the feature of interest
		record_focused = SeqRecord(
			id=record.id,
			annotations=record.annotations,
			dbxrefs=record.dbxrefs,
			seq=record.seq[start:end + 1],
			description=record.description
		)
		record_focused.features.append(highlight_feature)
	
		# Gather protein data for reference
		prep_prot_dict = {
		                  "nuccore_acc": record.id,
		                  # "region_seq": record.seq[start:end + 1],
		                  "window_start": start,
		                  "window_end": end,
		                  "feature_start": f_start,
		                  "feature_end": f_end,
		                  "strand": f_strand,
		                  "feature_len": f_len,
		                  "sequence": record.seq[start:end + 1]
		                  }
		prot_dict.setdefault(record.id,  prep_prot_dict)
		record_target.setdefault(record.id, record_focused)
	return prot_dict, record_target


def export_gbs(gb_dict, parent_path):
	if not os.path.exists(parent_path):
		os.mkdir(parent_path)
	for chromosome_id in gb_dict:
		query_suffix = re.sub(r'\|', '_', chromosome_id[0:20])
		out_path = f"{parent_path}{os.sep}{query_suffix}"
		if not os.path.exists(out_path):
			os.mkdir(out_path)

		gbk = gb_dict[chromosome_id]
		filename = f"{chromosome_id}.gb"
		with open(f"{out_path}{os.sep}{filename}", "w") as gb_handle:
			SeqIO.write(gbk, gb_handle, "genbank")


def main():
	# SNAKEMAKE IMPORTS
	# Inputs
	# Outputs
	# Params

	Entrez.email = config["entrez_login"]

	# Find NCBI's UID for a given HGVS query
	clinvar_uid = get_clinvar_uid('NM_001355224.2:c.867C>A')

	# Retrieve GENE object from Entrez based on clinvar
	(linked, hit_id, notfound) = elink_routine('clinvar', 'gene', clinvar_uid)

	# Generate a dictionary from the GENE object
	ncbi_gene_dict = get_gene_xml(linked)
	gene_coordinates = get_gene_coords(ncbi_gene_dict)

	(report_dict, focus_gb_dict) = get_genbank(gene_coordinates, 2000)

	export_gbs(focus_gb_dict, "jobs")


qual = []
for f in record.features:
	qual.append(f.qualifiers)

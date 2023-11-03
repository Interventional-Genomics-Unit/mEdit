# Native modules
import os
import copy
import pandas as pd
from pathlib import Path
# Installed modules
import urllib.error
import urllib3
import xmltodict
# Biopython
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Load config file
# with open("config/guide_prediction_default_template.yaml", "r") as f:
# 	config = yaml.load(f, Loader=yaml.FullLoader)


def list_lines_from_txt(file_path):
	with open(file_path, "r") as file_handle:
		list_of_entries = []
		for line in file_handle:
			clean_line = line.strip()
			list_of_entries.append(clean_line)
	return list_of_entries


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
	http = urllib3.PoolManager()
	file = http.request('GET', str(url))
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
	gene_uid = gene_xml_dict['Entrezgene-Set']['Entrezgene']['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']

	# Get gene coordinates for each individual chromosome
	for chr_index in range(0, chromosomes_list_length):
		chromosome_acc = chromosomes[chr_index]['Gene-commentary_accession']
		chromosome_version = chromosomes[chr_index]['Gene-commentary_version']
		chromosome_id = f"{chromosome_acc}.{chromosome_version}"
		gene_coords.setdefault(gene_uid, {}).setdefault(chromosome_id, [
			gene_xml_dict['Entrezgene-Set']['Entrezgene']['Entrezgene_locus']['Gene-commentary'][chr_index]['Gene-commentary_seqs']['Seq-loc']['Seq-loc_int']['Seq-interval']['Seq-interval_from'],
			gene_xml_dict['Entrezgene-Set']['Entrezgene']['Entrezgene_locus']['Gene-commentary'][chr_index]['Gene-commentary_seqs']['Seq-loc']['Seq-loc_int']['Seq-interval']['Seq-interval_to']
		])
	return gene_coords


def get_genbank(gene_coords, win_size):
	prot_dict = {}
	record_target = {}
	print("Fetch GenBank data via Entrez")
	for gene_uid in gene_coords:
		for chromosome_id in gene_coords[gene_uid]:
			print("Processing Gene gi", gene_uid, " from chromosome ", chromosome_id)
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
			f_start = int(gene_coords[gene_uid][chromosome_id][0])
			f_end = int(gene_coords[gene_uid][chromosome_id][1])
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

			# Create a fasta record with the feature of interes
			fasta_focused: record.seq[start:end + 1]

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
			                  }
			prot_dict.setdefault(gene_uid, {}).setdefault(record.id,  prep_prot_dict)
			record_target.setdefault(gene_uid, {}).setdefault(record.id, record_focused)
	return prot_dict, record_target


def create_new_dir(path):
	Path(path).mkdir(parents=True, exist_ok=True)


def export_gbs(gb_dict, hgvs_id, parent_path):
	for gene_uid in gb_dict:
		hgvs_output_directory = f"{parent_path}{os.sep}{gene_uid}:{hgvs_id}"
		genbank_output_directory = f"{hgvs_output_directory}{os.sep}gb"
		fasta_output_directory = f"{hgvs_output_directory}{os.sep}fna"
		create_new_dir(genbank_output_directory)
		create_new_dir(fasta_output_directory)
		for chromosome_id in gb_dict[gene_uid]:
			gbk = gb_dict[gene_uid][chromosome_id]
			with open(f"{genbank_output_directory}{os.sep}{chromosome_id}.gb", "w") as gb_handle:
				SeqIO.write(gbk, gb_handle, "genbank")
			with open(f"{fasta_output_directory}{os.sep}{chromosome_id}.fna", "w") as gb_handle:
				SeqIO.write(gbk, gb_handle, "fasta")


# DEBUG INPUT
# hgvs_id = 'NM_001355224.2:c.867C>A'
# entrez_login = 'thedoudnalab@gmail.com'
# genomic_window_size = 1000

def main():
	# SNAKEMAKE IMPORTS
	# Inputs
	hgvs_file_path = str(snakemake.input.hgvs_list)
	# Outputs
	report_out = str(snakemake.output.database_report)
	# Params
	genomic_window_size = int(snakemake.params.window_size)
	entrez_login = str(snakemake.params.entrez_login)
	# Wildcards
	output_directory = str(snakemake.wildcards.output_directory)
	run_name = str(snakemake.wildcards.run)

	# Get list of HGVSs from file
	hgvs_list = list_lines_from_txt(hgvs_file_path)

	# Set the analysis directory name
	analysis_directory = f"{output_directory}{os.sep}{run_name}"

	for hgvs_id in hgvs_list:
		# Entrez login
		Entrez.email = entrez_login

		# Find NCBI's UID for a given HGVS query
		clinvar_uid = get_clinvar_uid(hgvs_id)

		# Retrieve GENE object from Entrez based on clinvar
		(linked, hit_id, notfound) = elink_routine('clinvar', 'gene', clinvar_uid)

		# Generate a genomic coordinate dictionary from the GENE object
		ncbi_gene_dict = get_gene_xml(linked)
		gene_coordinates = get_gene_coords(ncbi_gene_dict)

		# Consolidate focused GenBank entries
		(report_dict, focus_gb_dict) = get_genbank(gene_coordinates, genomic_window_size)

		# Export focused GenBank to text files
		export_gbs(focus_gb_dict, hgvs_id, analysis_directory)

		# Export report file
		pd.DataFrame.from_dict(report_dict[list(report_dict.keys())[0]]).to_csv(report_out)


if __name__ == "__main__":
	main()

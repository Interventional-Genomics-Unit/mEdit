# Native Modules
import logging
import pickle
import subprocess
# Installed Modules
import pandas as pd
import vcf
import logging
# == Pysam REQUIRED ==
# Project Modules
from dataH import DataHandler

#################################
# Compares ALT VCF to Guide found in Hg38
# Searches and compares changes in either guide or pam
# If new guide is created --> submits query back to Fetch guides
#################################

def make_dict(ALTguides_dict, new_guides_set,new_beguides_set,set_header):

	if set_header:
		header = list(set(list(new_beguides_set.keys()) + list(new_guides_set.keys())))
		for k in header:
			ALTguides_dict[k] = []

	if len(new_guides_set['gRNA']) > 0:
		value_length = len(new_guides_set['gRNA'])
		for k in ALTguides_dict.keys():
			if k in new_guides_set.keys():
				ALTguides_dict[k] += [x for x in new_guides_set[k]]
			else:
				ALTguides_dict[k] += ["-"] * value_length

		if len(new_beguides_set['gRNA']) > 0:
			value_length = len(new_beguides_set['gRNA'])
			for k in ALTguides_dict.keys():
				if k in new_beguides_set.keys():
					ALTguides_dict[k] += [x for x in new_beguides_set[k]]
				else:
					ALTguides_dict[k] += ["-"] * value_length
		return ALTguides_dict


def extract_vcf_record(REFcoord, vcf_fname, window):
	'''
	searches for overlapping variant at the provided REF genome coordinate $ window
	'''
	records_found = []
	ch = REFcoord.split(':')[0]
	pos = REFcoord.split(':')[1].split('-')[0]
	vcf_reader = vcf.Reader(filename=vcf_fname, compressed=True)
	try:
		for record in vcf_reader.fetch(ch, int(pos) - window, int(pos) + window):
			records_found.append(record)
	except ValueError:
		pass
	return records_found


def make_ALTseq(extracted_seq, pos, ref, alt):
	'''
	Given the Reference Genome sequence (variant incorporated), adds mutations from ALT VCF to the sequence
	:param extracted_seq: Ref sequence
	:param pos: Alt genome variant start pos
	:param ref: Reference Allele
	:param alt: Alt Genome alt allele
	'''
	variant_seq = ''
	if len(ref) == len(alt):  ##substitution
		variant_seq = (extracted_seq[0:pos] + alt.lower() + extracted_seq[pos + len(alt):])
	else:
		logging.warning("SNV sites only accepted at this time. Allele {ref} will not be processed. Skipping")
	return variant_seq


def extract_variant_info(record, REFextracted_seq, REFcoord):
	'''
	for a given alt record found determine the ALT information
	and create a new extracted sequence incorporating the ALT allele
	'''
	## RECORD INFO FOR SAMPLE =1
	ALTalt_list, ALTref =record.ALT,record.REF
	call = record.samples[0]
	zyg = 'heterozygous' if call.is_het else 'homozygous'
	window = int(len(REFextracted_seq)/2)

	## if there is more than 2 ALT alleles
	if len(ALTalt_list) == 2:
		zyg = zyg + '-biallelic'
	if len(ALTalt_list) > 2:
		zyg = zyg + '-multiallelic'

	## THE first alt allele info
	vtype = record.ALT[0].type
	ALTalt = ALTalt_list[0].sequence ### <-- this program only incorporates the first alternative allele
	REFstart = int(REFcoord.split(':')[1]) - window
	ALTstart, ALTend = record.start, record.end
	ALTcoord = f'{record.CHROM}:{record.POS}'
	rel_pos = (ALTstart -REFstart)+ 1
	ALTalts = ",".join([str(x) for x in record.samples[0].site.alleles][1:])
	ALTid = f'{ALTcoord}:{ALTref}>{ALTalts}'

	#print('alt, ref',record.ALT, record.REF)
	#print('relative position',rel_pos)
	#print('check extracted seq == ref',REFextracted_seq[rel_pos],record.REF)
	#print('ALT starts and stop',ALTstart,ALTend)

	#if len(record.REF) == len(alt): ## substitution
	#    vtype = vtype + '-sub'
	#    new_seq = REFextracted_seq[0:rel_pos[0]] + alt.sequence + REFextracted_seq[rel_pos[1]:]
		#print(new_seq)

	#elif len(record.REF) > len(alt): # deletion
	#    vt = vtype + '-del'
		#new_seq = hg38extracted_seq[0:rel_pos[0]] + alt.sequence + hg38extracted_seq[rel_pos[1]:]
	#    new_seq = 'undetermined'

	#elif len(record.REF) < len(alt): # insertion
	#    vt = vtype + '-ins'
	#    new_seq = REFextracted_seq[0:rel_pos[0]] + alt.sequence + REFextracted_seq[rel_pos[1]:]
	#    new_seq = new_seq[0:len(REFextracted_seq)-1]


	#else:
	#    new_seq = 'undetermined'

	return ALTref, ALTalt, ALTalts, ALTcoord,zyg, vtype, rel_pos, ALTid


def find_overlapping_variants(filtered_vcf, models_path, hg38_snvinfo, search_params,be_search_params):
	ALTguides_dict = {}
	ALTinfo = {}
	set_header = True

	for ch, data in hg38_snvinfo.items():
		for d in data:
			query, tid, eid, gname,strand, REFref, REFalt, feature_annotation, REFextracted_seq, codons, REFcoord = d

			# Check VCF if variant exists in hg38 extracted_seq
			window = len(REFextracted_seq)/2
			records_found = extract_vcf_record(REFcoord=REFcoord, vcf_fname = filtered_vcf, window=window)

			if len(records_found) > 0:

				ALTseq = REFextracted_seq

				for rec in records_found:
					ALTref, ALTalt, ALTalts, ALTcoord, zyg, vtype, rel_pos, ALTid= extract_variant_info(rec,str(REFextracted_seq),REFcoord)

					if rel_pos != window: #make sure the alt genome mutation is not the same as the one in clinvar or it will be duplicated
						ALTinfo[query] = {ALTid : [ALTid,ALTref, ALTalt, ALTalts, ALTcoord, zyg, vtype, window -rel_pos if rel_pos>window else -1*rel_pos]}
						ALTseq = make_ALTseq(extracted_seq = ALTseq, pos = rel_pos, ref = REFalt, alt = ALTalt)

				if ALTseq != REFextracted_seq:
					# find quides
					dh = DataHandler(query, strand, REFref, REFalt, feature_annotation, models_path, ALTseq, codons,
									 REFcoord, gname)
					new_guides_set, new_beguides_set  = dh.get_Guides(search_params,be_search_params)

					make_dict(ALTguides_dict, new_guides_set, new_beguides_set, set_header)
					set_header = False
	return ALTinfo, ALTguides_dict


def guide_compare(guides_report, be_report, ALTguides_dict, ALTinfo, altgenome_name):
	'''
	Compares the guides generated from the ALT and REf genomes.
	Creates a DF that shows the ALT genomes impact on REF guides
	:param guides_report: Reference Genome guide CSV
	:param guides_report: Reference Genome base editor guide CSV
	:param ALTguides_dict: Guide found with variations found in ALT genome
	:param ALTinfo: ALT genome VCF variant into
	:param altgenome_name: example 'H02557'
	'''

	new_guides = []

	altgdf = pd.DataFrame(ALTguides_dict)
	refgdf = pd.read_csv(guides_report) # combine both be guides and endo guides together
	refbedf = pd.read_csv(be_report)
	refgdf = pd.concat([refgdf,refbedf])
	altgdf =altgdf.loc[:,refgdf.columns] # reorder cols to match

	refgdf = refgdf[refgdf['QueryTerm'].isin(list(ALTinfo.keys()))]
	REFrows = refgdf.to_dict('tight')['data']
	ALTrows = altgdf.to_dict('tight')['data']

	rename_columns = [f'Alt {x}' if 'core' in x or 'B' in x else x for x in altgdf.columns[8:]]
	columns = ['QueryTerm','GeneName', 'Editor', "Guide_ID",'Alt Guide Impact',"Alt Genome",
			   'Coordinates','Strand','Alt gRNA', 'Alt Pam',
			   f'Ref gRNA', f'Ref Pam'] + rename_columns

	for rrow in REFrows:
		query, gene, ed, gid, coord, strand, REFgrna, REFpam = rrow[0:8]

		cnt = 1
		for arow in ALTrows:
			if query == arow[0] and ed ==arow[2] :
				ALTgrna, ALTpam = arow[6:8]

				if [REFpam, REFgrna] == [ALTpam, ALTgrna]:  # unchanged
					ALTrows.remove(arow)
					cnt -= 1
					break

				elif REFpam == ALTpam and ALTgrna != REFgrna:  # same pam which means change is in grna
					newrow = rrow[0:3] + [f"alt_{rrow[3]}"]+ ['protospacer changed & conserved'] + [altgenome_name]+ rrow[4:6] + arow[6:8] + rrow[6:8] + arow[8:]
					new_guides.append(newrow)
					ALTrows.remove(arow)
					cnt -= 1
					break

				elif REFgrna == ALTgrna:  # pam is changed but conserved
					newrow = rrow[0:3] + [f"alt_{rrow[3]}"]+ ['PAM changed & conserved'] + [altgenome_name]+ rrow[4:6] + arow[6:8] + rrow[6:8] + arow[8:]
					new_guides.append(newrow)
					ALTrows.remove(arow)
					cnt -= 1
					break
				else:
					pass

		if cnt == 1:  # old guides remaining and not matched == no longer exsist
			newrow = rrow[0:3]  + [f"alt_{rrow[3]}"]+['PAM changed & removed'] +[altgenome_name]+  rrow[4:6] +["-","-"] + rrow[6:8] + len(rrow[8:]) * ['-']
			new_guides.append(newrow)


	if len(ALTrows) > 0:  # new guides remaining and not matched == new guides are made
		for arow in ALTrows:
			newrow = arow[0:3] + [f"alt_{arow[3]}_"] + ['PAM changed & added'] +[altgenome_name]+  arow[4:8] + ["-", "-"] + arow[8:]
			new_guides.append(newrow)

	ALTguides_df = pd.DataFrame(new_guides, columns=columns)
	return ALTguides_df


def create_ALTvcf_report(ALTinfo, altgenome_name):
	columns = ['QueryTerm', 'Genome', 'Alt VCF Variant ID',
			   'REF Allele', 'ALT Examined', 'ALT allele(s)',
			   'Variant Coordinates', 'Zygosity', 'Variant Type', 'Relative Position to Query']
	rows = []
	for query in ALTinfo.keys():
		# {ALTid : [ALTid,ALTref, ALTalt, ALTalts, ALTcoord, zyg, vtype, rel_pos]}
		for variant, info in ALTinfo[query].items():

			rows.append([query, altgenome_name] + info)

	nearby_variants_df = pd.DataFrame(rows, columns=columns)
	return nearby_variants_df


def fetch_ALT_guides(filtered_vcf,
					 guides_report,
					 be_report,
					 guide_search_params,
					 guide_be_search_params,
					 models_path,
					 snv_site_info,
					 diffguides_out,
					 altvar_out,
					 altgenome_name):

	# Get search parameters and the results from the reference assembly
	# {'spCas9': ('NGG', False, 20, -3, 'requirements work for SpCas9-HF1, eSpCas9 1.1,spyCas9'),
	search_params = pickle.load(open(guide_search_params, 'rb'))
	be_search_params = pickle.load(open(guide_be_search_params, 'rb'))
	if len(be_search_params.values()) == 0:
		be_search_params = None
	# get REF variant info
	# {'X': [['NM_004208.4:c.696+3G>A', 'NM_145812.3', '-', 'AIFM1', '-', 'C', 'T', 'intron', 'GGCTGG...
	hg38_snvinfo = pickle.load(open(snv_site_info, 'rb'))

	logging.info(f"Searching for changes in genome {altgenome_name}")
	# Check to see if there's alt variants in vcf and if so, find guides
	ALTinfo, ALTguides_dict = find_overlapping_variants(filtered_vcf, models_path, hg38_snvinfo, search_params,be_search_params)

	if len(ALTinfo.values()) > 0:

		nearby_variants_df = create_ALTvcf_report(ALTinfo,altgenome_name)
		nearby_variants_df.to_csv(altvar_out, index=False)
		if len(ALTguides_dict.values()) > 0:

			ALTguides_df = guide_compare(guides_report, be_report, ALTguides_dict, ALTinfo, altgenome_name)
			ALTguides_df.to_csv(diffguides_out, index=False)
		else:

			logging.info(f'No overlapping variants detected in {altgenome_name}')
			with open(diffguides_out, 'w') as f:
				f.write(f"No guide differences found based on the VCF {altgenome_name}")
	else:
		logging.info(f'no overlapping variants detected in {altgenome_name}')
		with open(altvar_out, 'w') as f:
			f.write(f"No Nearby variants found based on the VCF {altgenome_name}")


def main():
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	filtered_vcf = str(snakemake.input.filtered_vcf)
	guides_report = str(snakemake.input.guides_report_out)
	be_report = str(snakemake.input.be_report_out)
	guide_search_params = str(snakemake.input.guide_search_params)
	guide_be_search_params = str(snakemake.input.guide_be_search_params)
	snv_site_info = str(snakemake.input.snv_site_info)
	# === Outputs ===
	diffguides_out = str(snakemake.output.diff_guides)
	altvar_out = str(snakemake.output.alt_var)
	#   == Log File ==
	logfile_path = str(snakemake.output.logfile_path)
	# === Params ===
	models_path = str(snakemake.params.models_path)
	# === Wildcards ===
	altgenome_name = str(snakemake.wildcards.vcf_id)

	# === Log Process Initialization
	# Configure the logging system
	logging.basicConfig(
		level=logging.DEBUG,  # Set the minimum log level (DEBUG logs everything)
		format='%(asctime)s [%(levelname)s] %(message)s',  # Define log format
		handlers=[
			logging.FileHandler(logfile_path),  # Log to file
		]
	)

	logging.info('=== INITIALIZING ALTERNATIVE GENOME COMPARISON ROUTINE ===')

	# == Create dummy files for the outputs
	dummy_outputs = [diffguides_out, altvar_out]
	for report in dummy_outputs:
		# Create an empty DataFrame
		df = pd.DataFrame()
		# Export the empty DataFrame to a CSV file
		df.to_csv(str(report))

	fetch_ALT_guides(filtered_vcf, guides_report, be_report, guide_search_params, guide_be_search_params, models_path,
					 snv_site_info, diffguides_out, altvar_out, altgenome_name)

	logging.info('=== ALTERNATIVE GENOME COMPARISON ROUTINE FINALIZED ===')


if __name__ == "__main__":
	main()

# altgenome_name = "HG02257"
# outdir = "/groups/clinical/projects/editability/medit_queries/medit_test_3/standard/jobs/standard_TEST_standard/guide_prediction-hg38_GCA_000001405.15"
# db = "/groups/clinical/projects/editability/medit_queries/medit_database"
# filtered_vcf = f"/groups/clinical/projects/clinical_shared_data/hprc/hprc-v1.1-combined-phased-decomposed/{altgenome_name}.vcf.gz"
# guides_report = f"{outdir}/guides_report_ref/0_Guides_found.csv"
# be_report= f"{outdir}/guides_report_ref/0_BaseEditors_found.csv"
# guide_search_params = f"{outdir}/dynamic_params/0_guide_search_params.pkl"
# guide_be_search_params= f"{outdir}/dynamic_params/0_guide_be_search_params.pkl"
# models_path =f"{db}/pkl/models/"
# snv_site_info = f"{outdir}/dynamic_params/0_snv_site_info.pkl"
# diffguides_out = f"{outdir}/guides_report_{altgenome_name}/0_Guide_differences.csv"
# altvar_out = f"{outdir}/guides_report_{altgenome_name}/0_Alternative_genome_variants.csv"
# refgenome_name = "hg38"
#
#
# df =pd.read_csv(diffguides_out)
# df.Guide_ID

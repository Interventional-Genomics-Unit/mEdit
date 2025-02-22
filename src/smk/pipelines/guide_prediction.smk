# **** Import Packages ****
import glob

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# === Serialize Chromosomes for further processing
		expand("{root_dir}/tmp/{sequence_id}_chr_manifest.csv",
			root_dir=config["output_directory"],sequence_id=config["sequence_id"]),
		# === Predicted guides using the most recent human genome assembly ===
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_Guides_found.csv",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"],sequence_id=config["sequence_id"],
			query_index=config['query_index']),
		# === With the relevant VCF processed, proceed with creating consensus FASTA ===
		expand("{meditdb_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz",
			meditdb_path=config["meditdb_path"],mode=config["processing_mode"],
			vcf_id=config["vcf_id"],sequence_id=config["sequence_id"]),
		# === Predicted guides on alternative genomes based on the reference listed above ===
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/{query_index}_Guide_differences.csv",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"],sequence_id=config["sequence_id"],
			vcf_id=config["vcf_id"],query_index=config['query_index']),
		# === Compile Information from all alternative genomes (if any) into one single CSV ===
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_alt/{query_index}_Aggregated_Guide_differences.csv",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"],sequence_id=config["sequence_id"],
			query_index=config['query_index'])

# noinspection SmkAvoidTabWhitespace
rule serialize_chromosomes:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id))
	output:
		serialized_chr_manifest = "{root_dir}/tmp/{sequence_id}_chr_manifest.csv"
	params:
		decompressed_assembly=lambda wildcards: f"{config['fasta_root_path']}/{wildcards.sequence_id}.fa",
		output_dir=config["fasta_root_path"]
	conda:
		"../envs/tabix.yaml"
	threads:
		config["threads"]
	message:
		"""
# === SERIALIZE CHROMOSOMES ON REFERENCE GENOMES === #	
Inputs used:\n {input.assembly_path}
Outputs stored on:\n {params.output_dir}		
		"""
	script:
		"py/serialize_chromosomes.py"

# noinspection SmkAvoidTabWhitespace
rule predict_guides:
	input:
		query_manifest = "{root_dir}/queries/{run_name}_{query_index}.csv",
		serialized_chr_manifest = "{root_dir}/tmp/{sequence_id}_chr_manifest.csv",
		assembly_dir_path = lambda wildcards: glob.glob("{fasta_root_path}".format(
			fasta_root_path=config["fasta_root_path"]))
	output:
		guides_report_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_Guides_found.csv",
		nguides_report = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_Guides_count.csv",
		gene_report = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_Gene_Report.csv",
		variant_report = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_Variant_Report.csv",
		be_report_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_BaseEditors_found.csv",
		guide_search_params = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/{query_index}_guide_search_params.pkl",
		guide_be_search_params= "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/{query_index}_guide_be_search_params.pkl",
		snv_site_info = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/{query_index}_snv_site_info.pkl",
		guides_not_found_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_Guides_not-found.csv",
		logfile_path = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_guide_prediction.log",
	params:
		# == Main output path
		main_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref",
		# == Processed tables branch
		support_tables = config["support_tables"],
		annote_path = config["refseq_table"],
		# == Editor Parameters
		editors = config["editors"],
		base_editors = config["base_editors"],
		models_path= config["models_path"],
		distance_from_cutsite = config["distance_from_cutsite"],
		# == Custom Editor Parameters
		pam = config["pam"],
		guide_length = config["guide_length"],
		pam_is_first = config["pam_is_first"],
		dsb_position = config["dsb_position"],
		editing_window = config["editing_window"],
		target_base = config["target_base"],
		result_base = config["result_base"],
		# == Run Parameters ==
		qtype = config["qtype"],
		be_request = config["be_request"],
		editor_request = config["editor_request"],
	conda:
		"../envs/medit.yaml"
	message:
		"""
# === PREDICT GUIDES ON REFERENCE GENOMES === #	
Inputs used:
--> Take variants from:\n {input.query_manifest}
--> Use reference assembly:\n {input.assembly_dir_path}
--> Support tables from:\n {params.support_tables}

Run parameters:
--> Query type: {params.qtype} 
--> BEmode: {params.be_request}
--> Editor scope: {params.editor_request}

Outputs generated:
--> Generate reports on:\n {output.guides_report_out}\n {output.be_report_out}
Log file path:
--> {output.logfile_path}
        """
	script:
		"py/fetchGuides.py"

# noinspection SmkAvoidTabWhitespace
rule filter_vcf:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.filtered.fa.gz".format(
			fasta_root_path=config["fasta_root_path"], sequence_id=wildcards.sequence_id)),
		source_vcf=lambda wildcards: glob.glob("{meditdb_path}/{mode}/source_vcfs/{vcf_filename}.vcf.gz".format(
			meditdb_path=wildcards.meditdb_path, mode=wildcards.mode, vcf_filename=config["vcf_filename"]))
		# source_vcf="{meditdb_path}/{mode}/source_vcfs/{vcf_id}.vcf.gz"
	output:
		filtered_vcf="{meditdb_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz"
	params:
		split_vcf='{meditdb_path}/{mode}/source_vcfs/{vcf_id}.vcf.gz',
		source_vcf_prefix="{meditdb_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}",
		link_directory="{meditdb_path}/{mode}/consensus_refs/{sequence_id}/",
		dump_dir="{meditdb_path}/consensus_refs/downloads",
		fasta_root_path=config["fasta_root_path"]
	conda:
		"../envs/samtools.yaml"
	# resources:
	# 	mem_mb=100000
	message:
		"""
# === CREATE FILTERED VCF and CONSENSUS FASTA === #
This rule creates a consensus sequence based on a VCF file.
Inputs used:
--> Human genome assembly: {input.assembly_path}
--> Source VCF: {input.source_vcf}
Outputs generated:
--> Filtered VCF: {output.filtered_vcf}
Wildcards in this rule:
--> {wildcards}
		"""
	shell:
		"""
		# Prepare directories:
		# 1) If Depth is present in FORMAT filter > 5
		# 2) filter GAT1 or GAT2 samples (samples where one haplotype has a sequence depth = 0)
		# 3) filter reference & variant alleles > 1nt (SNVS only for now)
		# 4) Filter quality score > 15
		# Create index file
		
		# Inspect the number of genomes in the VCF file
		num_samples=$(bcftools query -l {input.source_vcf} | wc -l)

		# Check VCF only had one genome
		if [ $num_samples -le 1 ]; then
			# Check if read depth is given (for custom vcfs derived from wgs)
			if bcftools view -h {input.source_vcf} | grep -q '##FORMAT=<ID=DP'; then
				bcftools view --min-ac=1 -v snps -e 'GT="." || QUAL<15 || FMT/DP<5' {input.source_vcf} | bcftools sort -O z -o {output.filtered_vcf}
			else
				bcftools view --min-ac=1 -v snps -e 'GT="." || QUAL<15' {input.source_vcf} | bcftools sort -O z -o {output.filtered_vcf} 
			fi
			
		# If multi genomes split with sample name first
		else
			bcftools view -s {wildcards.vcf_id} --min-ac=1 {input.source_vcf} | bcftools sort -O z -o {output.filtered_vcf}
		fi

		bcftools index -f -t {output.filtered_vcf}		        		
        """

# noinspection SmkAvoidTabWhitespace
rule process_altgenomes:
	input:
		filtered_vcf=lambda wildcards: f"{config['meditdb_path']}/{wildcards.mode}/consensus_refs/{wildcards.sequence_id}/{wildcards.vcf_id}.filtered.vcf.gz",
		guides_report_out="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_Guides_found.csv",
		be_report_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_ref/{query_index}_BaseEditors_found.csv",
		guide_search_params="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/{query_index}_guide_search_params.pkl",
		guide_be_search_params= "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/{query_index}_guide_be_search_params.pkl",
		snv_site_info="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/{query_index}_snv_site_info.pkl"
	output:
		diff_guides = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/{query_index}_Guide_differences.csv",
		alt_var = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/{query_index}_Alternative_genome_variants.csv",
		logfile_path= "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/{query_index}_guide_prediction.log"
	params:
		idx_filtered_vcf = "{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz.tbi",
		models_path=config["models_path"],
	conda:
		"../envs/vcf.yaml"
	message:
		"""
# === PREDICT GUIDES ON ALTERNATIVE GENOMES === #	
Inputs used:
--> Template guides obtained from reference assembly:\n {input.guides_report_out}	
--> Processing guides based on VCF:\n {input.filtered_vcf}
--> Use reference assembly: {wildcards.sequence_id}
--> Take search parameters from:\n {input.guide_search_params}\n {input.snv_site_info}

Outputs generated:
--> Guide differences report output on:\n {output.diff_guides}
Wildcards in this rule:
--> {wildcards}
Log file path:
--> {output.logfile_path}
		"""
	script:
		"py/process_genome.py"

# noinspection SmkAvoidTabWhitespace
rule aggregate_altgenome_reports:
	input:
		diff_guides=expand("{{root_dir}}/{{mode}}/jobs/{{run_name}}/guide_prediction-{{sequence_id}}/guides_report_{vcf_id}/{{query_index}}_Guide_differences.csv",
			vcf_id=config["vcf_id"]),
		alt_var=expand("{{root_dir}}/{{mode}}/jobs/{{run_name}}/guide_prediction-{{sequence_id}}/guides_report_{vcf_id}/{{query_index}}_Alternative_genome_variants.csv",
			vcf_id=config["vcf_id"])
	output:
		aggr_diff_guides="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_alt/{query_index}_Aggregated_Guide_differences.csv",
		aggr_alt_var="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_alt/{query_index}_Alternative_genome_variants.csv",
		logfile_path= "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/guides_report_alt/{query_index}_guide_prediction.log"
	conda:
		"../envs/vcf.yaml"
	params:
		float_cols=config["float_cols"],
	message:
		"""
# === AGGREGATE REPORTS OF PREDICTED GUIDES ON ALTERNATIVE GENOMES === #	
Inputs used:
--> Guide differences report inputs:\n {input.diff_guides}
Outputs generated:
--> Aggregate Guide differences output:\n {output.aggr_diff_guides}
--> Aggregate Alternative variants output:\n {output.aggr_alt_var}
Log file path:
--> {output.logfile_path}
		"""
	script:
		"py/aggregate_alt_genomes.py"

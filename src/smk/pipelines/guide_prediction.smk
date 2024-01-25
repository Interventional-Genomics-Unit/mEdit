# **** Variables ****
# configfile: "config/guide_prediction_private_template.yaml"
# configfile: "config/preprocessing_configuration.yaml"

# configfile: "config/aws_download.yaml"

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile guide_prediction.smk -j 1 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --use-conda &

# Description:

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# With the relevant VCF downloaded, proceed with creating consensus FASTA
		expand("{meditdb_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.fa",
			meditdb_path=config["meditdb_path"],mode=config["processing_mode"],
			vcf_id=config["vcf_id"],sequence_id=config["sequence_id"]),
		# Predicted guides using the most recent human genome assembly
		expand("{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_ref/Guides_found.csv",
			root_dir=config["output_directory"], mode=config["processing_mode"],
			job_name=config["run_name"], sequence_id=config["sequence_id"]),
		# Predicted guides on alternative genomes based on the reference listed above
		expand("{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/Guide_differences.csv",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			job_name=config["run_name"],sequence_id=config["sequence_id"],
			vcf_id=config["vcf_id"])

# noinspection SmkAvoidTabWhitespace
rule consensus_fasta:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id)),
		source_vcf = "{meditdb_path}/{mode}/source_vcfs/{vcf_id}.vcf.gz"
	output:
		consensus_fasta = "{meditdb_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.fa",
		filtered_vcf = "{meditdb_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz"
	params:
		source_vcf_prefix="{meditdb_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}",
		dump_dir="{meditdb_path}/consensus_refs/downloads",
		fasta_root_path=config["fasta_root_path"]
	conda:
		"../envs/samtools.yaml"
	# resources:
	# 	mem_mb=100000
	message:
		"""
# === CREATE CONSENSUS FASTA === #
This rule creates a consensus sequence based on a VCF file.
Inputs used:
--> Human genome assembly: {input.assembly_path}
--> Source VCF: {input.source_vcf}
Outputs generated:
--> Filtered VCF: {output.filtered_vcf}
--> Consensus genome sequence: {output.consensus_fasta}
Wildcards in this rule:
--> {wildcards}
		"""
	shell:
		"""
		# Prepare directories:
        # 2) filter GAT1 or GAT2 samples (samples where one haplotype has a sequence depth = 0)
        # & filter reference & variant alleles > 5nt
        # Create index file
        bcftools filter -O z -o {output.filtered_vcf} -e 'GT="." || ILEN <= -5 || ILEN >= 5' {input.source_vcf} 
        bcftools index -f -t {output.filtered_vcf}

        # 3) Making a consensus
        #previously made a seperate hg38 Ref Fasta that only have standard chromsomes --> /groups/clinical/projects/editability/tables/raw_tables/VCFs/hg38_standard.fa.gz
        samtools dict {input.assembly_path} -o {params.fasta_root_path}/{wildcards.sequence_id}.dict
        samtools faidx {input.assembly_path} -o {input.assembly_path}.fai

        # gzip -dv {output.filtered_vcf}
        # bgzip {params.source_vcf_prefix}.filtered.vcf

        bcftools consensus -f {input.assembly_path} {output.filtered_vcf} -o {output.consensus_fasta}
        
        # Cleanup
        # rm {input.assembly_path}.fai {params.fasta_root_path}/{wildcards.sequence_id}.dict
        """

# noinspection SmkAvoidTabWhitespace
rule predict_guides:
	input:
		query_manifest = lambda wildcards: glob.glob("{variant_query_dir}".format(
			variant_query_dir=config["variant_query_dir"])),
		assembly_path = lambda wildcards: glob.glob("{fasta_root_path}".format(
			fasta_root_path=config["fasta_root_path"]))
	output:
		guides_report_out = "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_ref/Guides_found.csv",
		guide_search_params = "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/dynamic_params/guide_search_params.pkl",
		snv_site_info = "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/dynamic_params/snv_site_info.pkl"
	params:
		# == Main output path
		main_out = "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_ref",
		# == Main output filenames
		gene_report = config["gene_report"],
		variant_report = config["variant_report"],
		be_report = config["be_report"],
		# == Processed tables branch
		support_tables = config["support_tables"],
		annote_path = config["refseq_table"],
		# == Editor Parameters
		editors = config["editors"],
		base_editors = config["base_editors"],
		models_path= config["models_path"],
		# == Run Parameters ==
		qtype = config["qtype"],
		be_request = config["be_request"],
		editor_request = config["editor_request"]
	conda:
		"../envs/medit.yaml"
	message:
		"""
# === PREDICT GUIDES ON REFERENCE GENOMES === #	
Inputs used:
--> Take variants from:\n {input.query_manifest}
--> Use reference assembly:\n {input.assembly_path}
--> Support tables from:\n {params.support_tables}

Run parameters:
--> Query type: {params.qtype} 
--> BEmode: {params.be_request}
--> Editor scope: {params.editor_request}

Outputs generated:
--> Generate reports on: {output}
Wildcards in this rule:
--> {wildcards}
        """
	script:
		"py/fetchGuides.py"

# noinspection SmkAvoidTabWhitespace
rule process_altgenomes:
	input:
		filtered_vcf = lambda wildcards: glob.glob("{meditdb_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz".format(
			meditdb_path=config["meditdb_path"],mode=wildcards.mode,
			vcf_id=wildcards.vcf_id,sequence_id=wildcards.sequence_id
		)),
		guides_report_out= "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_ref/Guides_found.csv",
		guide_search_params= "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/dynamic_params/guide_search_params.pkl",
		snv_site_info= "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/dynamic_params/snv_site_info.pkl"
	output:
		diff_guides = "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/Guide_differences.csv",
	params:
		idx_filtered_vcf = "{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz.tbi",
		models_path= config["models_path"]
	# 	# == Main output path
	# 	main_out = "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/"
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
		"""
	script:
		"py/process_genome.py"

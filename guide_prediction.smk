# **** Variables ****
configfile: "config/guide_prediction_private_template.yaml"
configfile: "config/preprocessing_configuration.yaml"

# configfile: "config/aws_download.yaml"

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile guide_prediction.smk -j 1 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --use-conda &

# Description:

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Pull VCFs either from private (de novo sequenced) or the pangenomes available
		expand("{root_dir}/{mode}/source_vcfs/{vcf_id}.vcf.gz",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			vcf_id=config["vcf_id"]),
		# With the relevant VCF downloaded, proceed with creating consensus FASTA
		expand("{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.fa",
			root_dir=config["output_directory"],mode=config["processing_mode"],
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
rule pull_vcf:
	output:
		"{root_dir}/{mode}/source_vcfs/{vcf_id}.vcf.gz"
	params:
		aws_url = config["aws_url"],
		aws_path = config["aws_path"] ,
		aws_filename_suffix = config["filename_suffix"],
		vcf_id = config["vcf_id"],
		root_dir=config["output_directory"],
		mode=config["processing_mode"]
	shell:
		"""
        # 1) download diploid VCF files from AWS (-->to be a loop using index file) 
        touch {params.root_dir}/{params.mode}/source_vcfs/{wildcards.vcf_id}.vcf.gz
        wget {params.aws_url}/{wildcards.vcf_id}/{params.aws_path}/{wildcards.vcf_id}.{params.aws_filename_suffix}.vcf.gz -O {params.root_dir}/{params.mode}/source_vcfs/{wildcards.vcf_id}.vcf.gz || true
		"""

# noinspection SmkAvoidTabWhitespace
rule consensus_fasta:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id)),
		source_vcf = "{root_dir}/{mode}/source_vcfs/{vcf_id}.vcf.gz"
	output:
		consensus_fasta = "{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.fa",
		filtered_vcf = "{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz"
	params:
		source_vcf_prefix="{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}",
		dump_dir="{root_dir}/consensus_refs/downloads",
		fasta_root_path=config["fasta_root_path"]
	conda:
		"envs/samtools.yaml"
	resources:
		mem_mb=100000
	shell:
		"""
		# Prepare directories:
        # 2) filter GAT1 or GAT2 samples (samples where one haplotype has a sequence depth = 0)
        # & filter reference & variant alleles > 5nt
        # Create index file
        bcftools filter -O z -o {output.filtered_vcf} -e 'GT="." || ILEN <= -5 || ILEN >= 5' {input.source_vcf} 
        bcftools index -t {output.filtered_vcf}

        # 3) Making a consensus
        #previously made a seperate hg38 Ref Fasta that only have standard chromsomes --> /groups/clinical/projects/editability/tables/raw_tables/VCFs/hg38_standard.fa.gz
        samtools dict {input.assembly_path} -o {params.fasta_root_path}/{wildcards.sequence_id}.dict
        samtools faidx {input.assembly_path} -o {input.assembly_path}.fai

        gzip -dv {output.filtered_vcf}
        bgzip {params.source_vcf_prefix}.filtered.vcf

        bcftools consensus -f {input.assembly_path} {output.filtered_vcf} -o {output.consensus_fasta}

        # Cleanup
        rm {input.assembly_path}.fai {params.fasta_root_path}/{wildcards.sequence_id}.dict
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
		# == Run Parameters ==
		qtype = config["qtype"],
		BEmode = config["BEmode"],
		editor = config["editor"]
	conda:
		"envs/medit.yaml"
	message:
		"""
Take variants from:\n {input.query_manifest}
Run parameters:\n Query type: {params.qtype}; BEmode: {params.BEmode}; Editor scope: {params.editor}
Use reference assembly:\n {input.assembly_path}
Take support tables from:\n {params.support_tables}
Generate reports on:\n {output}
Wildcards: {wildcards}
        """
	script:
		"py/fetchGuides.py"

# noinspection SmkAvoidTabWhitespace
rule process_altgenomes:
	input:
		filtered_vcf = "{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz",
		guides_report_out= "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_ref/Guides_found.csv",
		guide_search_params= "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/dynamic_params/guide_search_params.pkl",
		snv_site_info= "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/dynamic_params/snv_site_info.pkl"
	output:
		diff_guides = "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/Guide_differences.csv",
	params:
		idx_filtered_vcf = "{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz.tbi"
	# 	# == Main output path
	# 	main_out = "{root_dir}/{mode}/jobs/{job_name}/guide_prediction-{sequence_id}/guides_report_{vcf_id}/"
	conda:
		"envs/vcf.yaml"
	message:
		"""
Template guides obtained from reference assembly:\n {input.guides_report_out}		
Processing guides based on VCF:\n {input.filtered_vcf}
Intermediate file generated by Tabix stored on:\n {params.idx_filtered_vcf}
Use reference assembly:\n {wildcards.sequence_id}
Take search parameters from:\n {input.guide_search_params}\n {input.snv_site_info}
Guide differences report output on:\n {output.diff_guides}
Wildcards: {wildcards}
		"""
	script:
		"py/process_genome.py"

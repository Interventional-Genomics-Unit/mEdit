# **** Variables ****
# configfile: "config/guide_prediction_default_template.yaml"
configfile: "config/guide_prediction_private_template.yaml"
# configfile: "config/aws_download.yaml"

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile vcf_processing.smk -j 1 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --use-conda &

# Description:

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Pull information from clinVar
		expand("{root_dir}/{mode}/source_vcfs/{vcf_id}.vcf.gz",
			root_dir=config["output_directory"], mode=config["processing_mode"], vcf_id=config["vcf_id"]),
		expand("{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.fa",
			root_dir=config["output_directory"], mode=config["processing_mode"],
			vcf_id=config["vcf_id"], sequence_id=config["sequence_id"])

rule pull_vcf:
	output:
		expand("{root_dir}/{mode}/source_vcfs/{vcf_id}.vcf.gz",
			root_dir=config["output_directory"], mode=config["processing_mode"], vcf_id=config["vcf_id"])
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
        wget {params.aws_url}/{params.vcf_id}/{params.aws_path}/{params.vcf_id}.{params.aws_filename_suffix}.vcf.gz -O {params.root_dir}/{params.mode}/source_vcfs/{params.vcf_id}.vcf.gz
		"""

rule consensus_fasta:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id)),
		source_vcf = "{root_dir}/{mode}/source_vcfs/{vcf_id}.vcf.gz"
	output:
		consensus_fasta = "{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.fa"
	params:
		source_vcf_prefix="{root_dir}/{mode}/source_vcfs/{vcf_id}",
		# target_consensus_prefix="{root_dir}/{mode}/consensus_refs/{sequence_id}/{vcf_id}",
		# genome_prefix=config["vcf_id"],
		dump_dir="{root_dir}/consensus_refs/downloads",
		fasta_root_path=config["fasta_root_path"]
	resources:
		mem_mb=100000
	shell:
		"""
		# Prepare directories:
        # 2) filter GAT1 or GAT2 samples (samples where one haplotype has a sequence depth = 0)
        # & filter reference & variant alleles > 5nt
        # Create index file
        bcftools filter -O z -o {params.source_vcf_prefix}.filtered.vcf.gz -e 'GT="." || ILEN <= -5 || ILEN >= 5' {input.source_vcf} 
        bcftools index -t {params.source_vcf_prefix}.filtered.vcf.gz

        # 3) Making a consensus
        #previously made a seperate hg38 Ref Fasta that only have standard chromsomes --> /groups/clinical/projects/editability/tables/raw_tables/VCFs/hg38_standard.fa.gz
        samtools dict {input.assembly_path} -o {params.fasta_root_path}/{wildcards.sequence_id}.dict
        samtools faidx {input.assembly_path} -o {input.assembly_path}.fai

        gzip -dv {params.source_vcf_prefix}.filtered.vcf.gz
        bgzip {params.source_vcf_prefix}.filtered.vcf

        bcftools consensus -f {input.assembly_path} {params.source_vcf_prefix}.filtered.vcf.gz -o {output.consensus_fasta}

        # Cleanup
        rm {input.assembly_path}.fai {input.assembly_path}.dict
        """
	# gatk IndexFeatureFile --input {output}/downloads/{params.genome_prefix}.filtered.vcf
	# gatk FastaAlternateReferenceMaker -R {input.assembly_path} -O {output}/fasta/{params.genome_prefix}.consensus.fa.gz -V {output}/downloads/{params.genome_prefix}.filtered.vcf

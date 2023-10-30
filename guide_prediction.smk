# **** Variables ****
configfile: "config/guide_prediction.yaml"
configfile: "config/aws_download.yaml"

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile guide_prediction.smk -j 1 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --use-conda &

# Description:

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Pull information from clinVar
		expand("{root_dir}/guide_prediction/{sequence_id}/guides_report/",
			sequence_id=config["sequence_id"],root_dir=config["output_directory"]),
		expand("{root_dir}/consensus_refs/{sequence_id}",
			sequence_id=config["sequence_id"],root_dir=config["output_directory"])

rule pull_vcf:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id))
	output:
		directory("{root_dir}/consensus_refs/{sequence_id}")
	params:
		genome_prefix=config["download_prefix"],
		dump_dir="{root_dir}/consensus_refs/downloads",
		fasta_root_path=config["fasta_root_path"]
	resources:
		mem_mb=100000
	shell:
		"""
        # Prepare directories:
        mkdir -p {output}/fasta {output}/downloads
        # 1) download diploid VCF files from AWS (-->to be a loop using index file) 
        wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/{params.genome_prefix}/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/{params.genome_prefix}.f1_assembly_v2.dip.vcf.gz -O {output}/downloads/{params.genome_prefix}.vcf.gz

        # 2) filter GAT1 or GAT2 samples (samples where one haplotype has a sequence depth = 0)
        # & filter reference & variant alleles > 5nt
        # Create index file
        bcftools filter -O z -o {output}/downloads/{params.genome_prefix}.filtered.vcf.gz -e 'GT="." || ILEN <= -5 || ILEN >= 5' {output}/downloads/{params.genome_prefix}.vcf.gz 
        bcftools index -t {output}/downloads/{params.genome_prefix}.filtered.vcf.gz

        # 3) Making a consensus
        #previously made a seperate hg38 Ref Fasta that only have standard chromsomes --> /groups/clinical/projects/editability/tables/raw_tables/VCFs/hg38_standard.fa.gz
        samtools dict {input.assembly_path} -o {params.fasta_root_path}/{wildcards.sequence_id}.dict
        samtools faidx {input.assembly_path} -o {input.assembly_path}.fai

        gzip -dv {output}/downloads/{params.genome_prefix}.filtered.vcf.gz
        bgzip {output}/downloads/{params.genome_prefix}.filtered.vcf
        cp {output}/downloads/{params.genome_prefix}.filtered.vcf.gz .

        bcftools consensus -f {input.assembly_path} {output}/downloads/{params.genome_prefix}.filtered.vcf.gz -o {output}/fasta/{params.genome_prefix}.consensus.fa
        # gatk IndexFeatureFile --input {output}/downloads/{params.genome_prefix}.filtered.vcf    
        # gatk FastaAlternateReferenceMaker -R {input.assembly_path} -O {output}/fasta/{params.genome_prefix}.consensus.fa.gz -V {output}/downloads/{params.genome_prefix}.filtered.vcf

        # Cleanup
        rm {input.assembly_path}.fai {input.assembly_path}.dict
        """


rule fetch_guides:
	#
	input:
		query_manifest=lambda wildcards: glob.glob("{variant_query_dir}/hgvs_test_queries.csv".format(
			variant_query_dir=config["variant_query_dir"])),
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id))
	output:
		directory("{root_dir}/guide_prediction/{sequence_id}/guides_report")
	params:
		support_tables=config["support_tables"]
	conda:
		"envs/medit.yaml"
	message:
		"""
        Take variants from:\n {input.query_manifest}
        Use reference assembly:\n {input.assembly_path}
        Take support tables from:\n {params.support_tables}
        Generate reports on:\n {output}
        Wildcards: {wildcards}
        """
	script:
		"py/fetchGuides.py"

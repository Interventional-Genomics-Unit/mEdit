# **** Variables ****
# configfile: "config/guide_prediction_default_template.yaml"
# configfile: "config/guide_prediction_private_template.yaml"
# configfile: "config/aws_download.yaml"

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile vcf_processing.smk -j 1 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --use-conda &

# Description:

# noinspection SmkAvoidTabWhitespace
# rule all:
# 	input:
# 		# Pull VCFs either from private (de novo sequenced) or the pangenomes available
# 		expand("{db_path}/{mode}/source_vcfs/{vcf_id}.vcf.gz",
# 			db_path=config["db_path"],mode=config["processing_mode"],
# 			vcf_id=config["vcf_id"]),
# 		# With the relevant VCF downloaded, proceed with creating consensus FASTA
# 		expand("{db_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.fa",
# 			db_path=config["db_path"],mode=config["processing_mode"],
# 			vcf_id=config["vcf_id"],sequence_id=config["sequence_id"])
#
# # noinspection SmkAvoidTabWhitespace
# rule pull_vcf:
# 	output:
# 		"{db_path}/{mode}/source_vcfs/{vcf_id}.vcf.gz"
# 	params:
# 		aws_url = config["aws_url"],
# 		aws_path = config["aws_path"] ,
# 		aws_filename_suffix = config["filename_suffix"],
# 		vcf_id = config["vcf_id"],
# 		db_path=config["db_path"],
# 		mode=config["processing_mode"]
# 	shell:
# 		"""
#         # 1) download diploid VCF files from AWS (-->to be a loop using index file)
#         touch {params.db_path}/{params.mode}/source_vcfs/{wildcards.vcf_id}.vcf.gz
#         wget {params.aws_url}/{wildcards.vcf_id}/{params.aws_path}/{wildcards.vcf_id}.{params.aws_filename_suffix}.vcf.gz -O {params.db_path}/{params.mode}/source_vcfs/{wildcards.vcf_id}.vcf.gz || true
# 		"""
#
# # noinspection SmkAvoidTabWhitespace
# rule consensus_fasta:
# 	input:
# 		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
# 			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id)),
# 		source_vcf = "{db_path}/{mode}/source_vcfs/{vcf_id}.vcf.gz"
# 	output:
# 		consensus_fasta = "{db_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.fa",
# 		filtered_vcf = "{db_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}.filtered.vcf.gz"
# 	params:
# 		source_vcf_prefix="{db_path}/{mode}/consensus_refs/{sequence_id}/{vcf_id}",
# 		dump_dir="{db_path}/consensus_refs/downloads",
# 		fasta_root_path=config["fasta_root_path"]
# 	conda:
# 		"envs/samtools.yaml"
# 	# resources:
# 	# 	mem_mb=100000
# 	shell:
# 		"""
# 		# Prepare directories:
#         # 2) filter GAT1 or GAT2 samples (samples where one haplotype has a sequence depth = 0)
#         # & filter reference & variant alleles > 5nt
#         # Create index file
#         bcftools filter -O z -o {output.filtered_vcf} -e 'GT="." || ILEN <= -5 || ILEN >= 5' {input.source_vcf}
#         bcftools index -t {output.filtered_vcf}
#
#         # 3) Making a consensus
#         #previously made a seperate hg38 Ref Fasta that only have standard chromsomes --> /groups/clinical/projects/editability/tables/raw_tables/VCFs/hg38_standard.fa.gz
#         samtools dict {input.assembly_path} -o {params.fasta_root_path}/{wildcards.sequence_id}.dict
#         samtools faidx {input.assembly_path} -o {input.assembly_path}.fai
#
#         gzip -dv {output.filtered_vcf}
#         bgzip {params.source_vcf_prefix}.filtered.vcf
#
#         bcftools consensus -f {input.assembly_path} {output.filtered_vcf} -o {output.consensus_fasta}
#         # --> DEBUG <--
# 		touch {output.consensus_fasta}
# 		# --> DEBUG <--
#         # Cleanup
#         rm {input.assembly_path}.fai {params.fasta_root_path}/{wildcards.sequence_id}.dict
#         """

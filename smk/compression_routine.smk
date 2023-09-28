# **** Variables ****
import glob

configfile: "../config/guide_prediction.yaml"

# Cluster run template
# nohup snakemake --snakefile compression_routine.smk -j 10 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config ../config/cluster.yaml --use-conda --use-conda &

# Description:


# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        expand("{fasta_download_path}/{assembly_id}/{sequence_id}.maternal.f1_assembly_v2_genbank.fa",
            assembly_id=config["assembly_id"],sequence_id=config["sequence_id"],fasta_download_path=config["fasta_download_path"]),
        expand("{fasta_download_path}/{assembly_id}/{sequence_id}.maternal.f1_assembly_v2_genbank.fa.gz",
               assembly_id=config["assembly_id"], sequence_id=config["sequence_id"], fasta_download_path=config["fasta_download_path"])

# rule aws_download:
#     input:
#         seq_id = lambda wildcards: glob.glob(wildcards.sequence_id)
#     output:
#         mat_aws = "{fasta_download_path}/{assembly_id}/aws/{sequence_id}.maternal.f1_assembly_v2_genbank.fa.gz",
#         pat_aws = "{fasta_download_path}/{assembly_id}/aws/{sequence_id}.paternal.f1_assembly_v2_genbank.fa.gz"
#     params:
#         aws_path = config["aws_s3_path"],
#         fasta_download_path = config["fasta_download_path"]
#     shell:
#         """
#         aws s3 cp {params.aws_path}/{wildcards.sequence_id}.maternal.f1_assembly_v2_genbank.fa.gz {params.fasta_download_path}/{wildcards.assembly_id}/aws/
#         aws s3 cp {params.aws_path}/{wildcards.sequence_id}.paternal.f1_assembly_v2_genbank.fa.gz {params.fasta_download_path}/{wildcards.assembly_id}/aws/
#         """


rule pigz_decompress:
    input:
        mat_aws="{fasta_download_path}/aws/{sequence_id}.maternal.f1_assembly_v2_genbank.fa.gz",
        pat_aws="{fasta_download_path}/aws/{sequence_id}.paternal.f1_assembly_v2_genbank.fa.gz"
    output:
        dcmp_mat = "{fasta_download_path}/{sequence_id}.maternal.f1_assembly_v2_genbank.fa",
        dcmp_pat= "{fasta_download_path}/{sequence_id}.paternal.f1_assembly_v2_genbank.fa"
    threads:
        config["threads"]
    shell:
        """
        pigz -dv -p {threads} {input.mat_aws}
        pigz -dv -p {threads} {input.pat_aws}
        """
    
rule bgzip_compress:
    input:
        dcmp_mat = "{fasta_download_path}/{sequence_id}.maternal.f1_assembly_v2_genbank.fa",
        dcmp_pat = "{fasta_download_path}/{sequence_id}.paternal.f1_assembly_v2_genbank.fa"
    output:
        bgz_mat = "{fasta_download_path}/{sequence_id}.maternal.f1_assembly_v2_genbank.fa.gz",
        bgz_pat = "{fasta_download_path}/{sequence_id}.paternal.f1_assembly_v2_genbank.fa.gz"
    threads:
        config["threads"]
    shell:
        """
        bgzip --threads {threads} {input.dcmp_mat}
        bgzip --threads {threads} {input.dcmp_pat}
        """

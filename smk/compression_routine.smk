# **** Variables ****
import glob

configfile: "../config/guide_prediction.yaml"
configfile: "../config/aws_download.yaml"

# Cluster run template
# nohup snakemake --snakefile compression_routine.smk -j 10 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config ../config/cluster.yaml --use-conda --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        expand("{fasta_compressed_path}/{sequence_id}{filename_suffix}.fa",
            fasta_compressed_path=config["fasta_compressed_path"],sequence_id=config["sequence_id"],filename_suffix=config["filename_suffix"]),
        expand("{fasta_compressed_path}/{sequence_id}{filename_suffix}.fa.gz",
            fasta_compressed_path=config["fasta_compressed_path"],sequence_id=config["sequence_id"],filename_suffix=config["filename_suffix"])

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
        downloaded_file = lambda wildcards: glob.glob("{fasta_download_path}/{sequence_id}{filename_suffix}.fa.gz".format(
            fasta_download_path=config["fasta_download_path"],sequence_id=wildcards.sequence_id,filename_suffix=wildcards.filename_suffix
        ))
    output:
        decompressed = "{fasta_compressed_path}/{sequence_id}{filename_suffix}.fa"
    params:
        download_path = config["fasta_download_path"]
    threads:
        config["threads"]
    shell:
        """
        pigz -dvk -p {threads} {input.downloaded_file}
        mv {params.download_path}/{wildcards.sequence_id}{wildcards.filename_suffix}.fa {wildcards.fasta_compressed_path}/{wildcards.sequence_id}{wildcards.filename_suffix}.fa
        sleep 1m
        """
    
rule bgzip_compress:
    input:
        decompressed = "{fasta_compressed_path}/{sequence_id}{filename_suffix}.fa"
    output:
        bgz_compressed = "{fasta_compressed_path}/{sequence_id}{filename_suffix}.fa.gz"
    threads:
        config["threads"]
    shell:
        """
        bgzip {input.decompressed}
        """

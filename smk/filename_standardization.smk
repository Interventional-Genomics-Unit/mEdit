# **** Variables ****
import glob

configfile: "../config/aws_download.yaml"
configfile: "../config/guide_prediction.yaml"

# Cluster run template
# nohup snakemake --snakefile filename_standardization.smk -j 10 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config ../config/cluster.yaml --use-conda --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        expand("{fasta_root_path}/{sequence_id}.fa.gz",
            sequence_id=config["sequence_id"],fasta_root_path=config["fasta_root_path"]),

rule symbolic_link:
    input:
        filename = lambda wildcards: glob.glob("{fasta_compressed_path}/{sequence_id}{filename_suffix}.fa.gz".format(
            fasta_compressed_path=config["fasta_compressed_path"],sequence_id=wildcards.sequence_id,filename_suffix=config["filename_suffix"]
        ))
    output:
        symlink_name = "{fasta_root_path}/{sequence_id}.fa.gz"
    shell:
        """
        ln -s {input.filename} {output.symlink_name}
        """

# **** Variables ****
import glob

configfile: "../config/guide_prediction_default_template.yaml"

# Cluster run template
# nohup snakemake --snakefile filename_standardization.smk -j 10 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config ../config/cluster.yaml --use-conda --latency-wait 120 &

# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        expand("{fasta_root_path}/{sequence_id}.fa.gz",
            fasta_root_path=config["fasta_root_path"], sequence_id=config["sequence_id"]),

rule symbolic_link:
    input:
        filename = lambda wildcards: glob.glob("{fasta_compressed_path}/{sequence_id}.{filename_suffix}.fa.gz".format(
            fasta_compressed_path=config["fasta_compressed_path"],
            sequence_id=wildcards.sequence_id,
            filename_suffix=config["filename_suffix"]))
    output:
        symlink_name = "{fasta_root_path}/{sequence_id}.fa.gz"
    params:
        symlink = "{sequence_id}.fa.gz"
    message:
        """
        Convert:\n {input.filename}
        To:\n {output.symlink_name}
        Available Wildcards: {wildcards}
        """
    shell:
        """
        cd {wildcards.fasta_root_path}
        ln -s {input.filename} {params.symlink}
        """

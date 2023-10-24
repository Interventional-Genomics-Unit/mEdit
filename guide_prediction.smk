# **** Variables ****
configfile: "config/guide_prediction.yaml"
# configfile: "config/hgvs_input.yaml"

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile guide_prediction.smk -j 1 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --use-conda --latency-wait 120 &

# Description:

# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        # Pull information from clinVar
        expand("{root_dir}/{sequence_id}/guides_report/",
            sequence_id=config["sequence_id"],root_dir=config["output_directory"])


rule select_assembly:



rule fetch_guides:
    #
    input:
        query_manifest = lambda wildcards: glob.glob("{variant_query_dir}/hgvs_test_queries.csv".format(
            variant_query_dir=config["variant_query_dir"])),
        assembly_path = lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
            fasta_root_path=config["fasta_root_path"], sequence_id=wildcards.sequence_id)),
        annote_path = lambda wildcards: glob.glob("{support_tables}/processed_tables/ncbiRefSeq.txt.gz".format(
            support_tables=config["support_tables"]))
    output:
        directory("{root_dir}/{sequence_id}/guides_report")
    params:
        support_tables = config["support_tables"]
    conda:
        "envs/medit.yaml"
    message:
        """
        Take variants from:\n {input.query_manifest}
        Use reference assembly:\n {input.assembly_path}
        Use support annotation:\n {input.annote_path}
        Take support tables from:\n {params.support_tables}
        Generate reports on:\n {output}
        Wildcards: {wildcards}
        """
    script:
        "py/fetchGuides.py"

# **** Variables ****
configfile: "config/guide_prediction.yaml"
configfile: "config/hgvs_input.yaml"

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile guide_prediction.smk -j 1 --cluster "sbatch -t {cluster.time}
# -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --use-conda &

# Description:


# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        # Pull information from clinVar
        expand("{root_dir}/{assembly_id}/guides_report/",
            assembly_id=config["assembly_id"],root_dir=config["root_dir"])

rule fetch_guides:
    #
    input:
        query_manifest = "{root_dir}/test_in/hgvs_test_queries.csv",
        assembly_path = lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
            fasta_root_path=config["fasta_root_path"], sequence_id=config["sequence_id"]))
    output:
        directory("{root_dir}/{sequence_id}/guides_report/")
    params:
        support_tables = config[""]
    conda:
        ""
    script:
        "../py/fetchGuides.py"

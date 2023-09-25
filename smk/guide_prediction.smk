# **** Variables ****
configfile: "config/guide_prediction.yaml"
configfile: "config/hgvs_input.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile guide_prediction.smk -j 1 --cluster "sbatch -t {cluster.time}
# -n {cluster.cores} -N {cluster.nodes}" --cluster-config config/cluster.yaml --use-conda &

# Description:


# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        # Pull information from clinVar
        expand("{output_directory}/{run}/clinvar_pull.report",
            output_directory=config["output_directory"],
            root_dir=config["root_dir"])

rule fetch_guides:
    #
    input:
        query_manifest = "{root_dir}/test_in/hgvs_test_queries.csv",
        assembly_path = "{root_dir}/{assembly_id}"
    output:
        directory("{root_dir}/guides_report/")
    params:
        window_size = config["window_size"],
        entrez_login = config["entrez_login"]
    conda:
        "../envs/clinvar_pull.yaml"
    script:
        "../py/fetchGuides.py"

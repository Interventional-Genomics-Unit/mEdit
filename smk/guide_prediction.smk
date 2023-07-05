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
            run=config["run_name_suffix"])

rule clinVar_pull:
    # Pulls genomic context, based on the {window_size}, of a gene(s) associated with a
    #   HGVS (or a list of) identifier(s).
    input:
        hgvs_list = config["hgvs_path"]
    output:
        database_report = "{output_directory}/{run}/clinvar_pull.report"
    params:
        window_size = config["window_size"],
        entrez_login = config["entrez_login"]
    conda:
        "../envs/clinvar_pull.yaml"
    script:
        "../py/clinVar_pull.py"

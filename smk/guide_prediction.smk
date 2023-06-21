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
        expand("{output_directory}/{hgvs}/clinvar_pull.report",
            output_directory=config["output_directory"],
            hgvs=config["hgvs_list"])

rule clinVar_pull:
    input:
        # TODO: Aqui vai precisar de um arquivo texto contendo os HGVSs
        hgvs_list = config["hgvs_list"]
    output:
        # TODO: Isso feito, a estrutura de nomes do output vai ter q mudar tb
        database_report = "{output_directory}/{hgvs}/clinvar_pull.report"
    params:
        window_size = config["window_size"],
        entrez_login = config["entrez_login"]
    conda:
        ""
    script:
        "../py/clinVar_pull.py"

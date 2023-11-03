# **** Variables ****
configfile: "config/guide_prediction_default_template.yaml"
configfile: "config/preprocessing_configuration.yaml"

# configfile: "config/aws_download.yaml"

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile guide_prediction.smk -j 1 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --use-conda &

# Description:

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Pull information from clinVar
		expand("{root_dir}/{mode}/jobs/{job_name}/guide_prediction/{sequence_id}/guides_report/",
			root_dir=config["output_directory"], mode=config["processing_mode"],
			job_name=config["run_name"], sequence_id=config["sequence_id"])


rule fetch_guides:
	input:
		query_manifest=lambda wildcards: glob.glob("{variant_query_dir}/hgvs_test_queries.csv".format(
			variant_query_dir=config["variant_query_dir"])),
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id))
	output:
		directory("{root_dir}/{mode}/jobs/{job_name}/guide_prediction/{sequence_id}/guides_report")
	params:
		support_tables=config["support_tables"]
	conda:
		"envs/medit.yaml"
	message:
		"""
        Take variants from:\n {input.query_manifest}
        Use reference assembly:\n {input.assembly_path}
        Take support tables from:\n {params.support_tables}
        Generate reports on:\n {output}
        Wildcards: {wildcards}
        """
	script:
		"py/fetchGuides.py"

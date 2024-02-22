# **** Variables ****
configfile: ""

# **** Imports ****
import glob

# Cluster run template
# nohup snakemake --snakefile *.smk -j 1 --cluster "sbatch -t {cluster.time} -n {cluster.cores}" --cluster-config config/cluster.yaml --use-conda &

# Description:

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Description
		expand("",),
		# Description
		expand("",),

# noinspection SmkAvoidTabWhitespace
rule casoff_input_formatting:
	input:
		guides_per_editor_path = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}.pkl",
		guide_search_params = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/guide_search_params.pkl",
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id))
	output:
		casoff_input = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff_in.txt",
		seq_pam_path = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_seqpam.pkl"
	params:
		tmp_processing_casoff = "",
		RNAbb = config["RNAbb"],
		DNAbb= config["DNAbb"],
		mm= config["mm"],
		PU = config["PU"]
	conda:
		"envs/"
	message:
		"""
		"""
	script:
		"py/build_casoff_input.py"

# noinspection SmkAvoidTabWhitespace
rule casoff_run:
	input:
		guides_report_out = ""
	output:
		casoff_out = ""
	params:
		tmp_casoff = "",
		RNAbb=config["RNAbb"],
		DNAbb=config["DNAbb"],
		mm=config["mm"],
		PU=config["PU"]
	conda:
		"envs/"
	threads:
		config["threads"]
	message:
		"""
		"""
	shell:
		"""
		"""

rule casoff_output_formatting:
	input:
		guides_report_out=""
	output:
		casoff_out = ""
	params:
		tmp_casoff = ""
	conda:
		"envs/"
	message:
		"""
		"""
	script:
		"py/"


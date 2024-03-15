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
		expand("{fasta_root_path}/{sequence_id}.fa",
			fasta_root_path=config["fasta_root_path"], sequence_id=config["sequence_id"]),
		# Prepare input files for casoffinder on a per-editor basis
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff_in.txt",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"], sequence_id=config["sequence_id"],
			editing_tool=config["editors_list"]),


rule decompress_genome:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id)),
	output:
		decompressed_assembly_symlink = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{sequence_id}.fa",
	params:
		decompressed_assembly_path = lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id)),
		link_directory = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/"
	priority: 50
	message:
		"""
		"""
	shell:
		"""
		gzip -kdv {input.assembly_path}
		ln --symbolic -t {params.link_directory} {params.decompressed_assembly_path}
		"""

# noinspection SmkAvoidTabWhitespace
rule casoff_input_formatting:
	input:
		guides_per_editor_path = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}.pkl",
		guide_search_params = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/guide_search_params.pkl",
		decompressed_assembly_symlink = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{sequence_id}.fa"
	output:
		casoff_input = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff_in.txt",
		seq_pam_path = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_seqpam.pkl"
	params:
		tmp_processing_casoff = config["tmp_processing_casoff"],
		rna_bulge = config["RNAbb"],
		dna_bulge= config["DNAbb"],
		max_mismatch= config["max_mismatch"],
		casoff_accelerator = config["PU"]
	conda:
		"../envs/casoff.yaml"
	message:
		"""
# === PREDICT OFFTARGET EFFECT === #	
Inputs used:
--> Take guides grouped by editing tool:\n {input.guides_per_editor_path}
--> Use reference assembly:\n {input.decompressed_assembly_symlink}
--> Use guide search parameters from:\n {input.guide_search_params}

Run parameters:
--> RNA bulge: {params.rna_bulge} 
--> DNA bulge: {params.dna_bulge}
--> Maximum mismatch: {params.max_mismatch}

Outputs generated:
--> CasOffinder formatted input: {output.casoff_input}
Wildcards in this rule:
--> {wildcards}
		"""
	script:
		"py/build_casoff_input.py"

# noinspection SmkAvoidTabWhitespace
rule casoff_run:
	input:
		casoff_input="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff_in.txt",
		seq_pam_path="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_seqpam.pkl",
		guides_report_out = ""
	output:
		casoff_out = ""
	params:
		tmp_processing_casoff=config["tmp_processing_casoff"],
		rna_bulge=config["RNAbb"],
		dna_bulge=config["DNAbb"],
		max_mismatch=config["max_mismatch"],
		casoff_accelerator=config["PU"]
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

# rule casoff_output_formatting:
# 	input:
# 		guides_report_out=""
# 	output:
# 		casoff_out = ""
# 	params:
# 		tmp_casoff = ""
# 	conda:
# 		"envs/"
# 	message:
# 		"""
# 		"""
# 	script:
# 		"py/"


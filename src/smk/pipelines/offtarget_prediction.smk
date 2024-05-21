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
		# expand("{fasta_root_path}/{sequence_id}.fa",
		# 	fasta_root_path=config["fasta_root_path"], sequence_id=config["sequence_id"]),
		# Prepare input files for casoffinder on a per-editor basis
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/input_files/{editing_tool}_casoff_in.txt",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"], sequence_id=config["sequence_id"],
			editing_tool=config["editors_list"]),
		# Run Cas-Offinder
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff.txt",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"], sequence_id=config["sequence_id"],
			editing_tool=config["editors_list"]),
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff_parsed.txt",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"], sequence_id=config["sequence_id"],
			editing_tool=config["editors_list"])


rule decompress_genome:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id))
	output:
		decompressed_assembly_symlink = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{sequence_id}.fa"
	params:
		decompressed_assembly_path = lambda wildcards: glob.glob("{fasta_root_path}/{sequence_id}.fa".format(
			fasta_root_path=config["fasta_root_path"],sequence_id=wildcards.sequence_id)),
		link_directory = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/"
	priority: 50
	message:
		"""
		# === PREPARING REFERENCE GENOMES FOR CASOFFINDER === #
		Inputs used:
		Compressed genome reference: {input.assembly_path}
		Outputs:
		Target directoru: {params.link_directory}
		Decompressed genome reference: {output.decompressed_assembly_symlink}
		Wildcards in this rule:
		{wildcards}
		"""
	shell:
		"""
		gzip -kdvf {input.assembly_path}
		ln --symbolic -t {params.link_directory} {params.decompressed_assembly_path}
		"""

# noinspection SmkAvoidTabWhitespace
rule casoff_input_formatting:
	input:
		guides_per_editor_path = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}.pkl",
		guide_search_params = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/dynamic_params/guide_search_params.pkl",
		decompressed_assembly_symlink = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{sequence_id}.fa",
	output:
		casoff_input = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/input_files/{editing_tool}_casoff_in.txt",
		casoff_support = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/input_files/{editing_tool}_casoff_support.txt"
	params:
		rna_bulge = config["RNAbb"],
		dna_bulge= config["DNAbb"],
		max_mismatch= config["max_mismatch"],
		casoff_accelerator = config["PU"]
	conda:
		"../envs/casoff.yaml"
	message:
		"""
# === DATA FORMATTING FOR CAS-OFFINDER === #	
Inputs used:
--> Take guides grouped by editing tool:\n {input.guides_per_editor_path}
--> Use reference assembly:\n {input.decompressed_assembly_symlink}
--> Use guide search parameters from:\n {input.guide_search_params}
--> Temp files stored at:\n {output.casoff_support}

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
		casoff_input="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/input_files/{editing_tool}_casoff_in.txt"
	output:
		casoff_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff.txt"
	params:
		rna_bulge=config["RNAbb"],
		dna_bulge=config["DNAbb"],
		max_mismatch=config["max_mismatch"],
		casoff_accelerator=config["PU"]
	conda:
		"../envs/casoff.yaml"
	threads:
		int(config["threads"])
	message:
		"""
# === PREDICT OFFTARGET EFFECT === #
Inputs used:
--> Analyze off-target effect for guides predicted for: {wildcards.editing_tool}
--> Take formatted inputs from :\n {input.casoff_input}

Run parameters:
--> RNA bulge: {params.rna_bulge} 
--> DNA bulge: {params.dna_bulge}
--> Maximum mismatch: {params.max_mismatch}
--> Cas-Offinder running on device: {params.casoff_accelerator}

Outputs generated:
--> CasOffinder output: {output.casoff_out}
Wildcards in this rule:
--> {wildcards}		
		"""
	shell:
		"""
		cas-offinder {input.casoff_input} {params.casoff_accelerator} {output.casoff_out}
		"""

# noinspection SmkAvoidTabWhitespace
rule casoff_output_formatting:
	input:
		casoff_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff.txt",
		casoff_support = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/input_files/{editing_tool}_casoff_support.txt",
	output:
		formatted_casoff_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{sequence_id}/offtarget_prediction/{editing_tool}_casoff_parsed.txt"
	params:
		rna_bulge=config["RNAbb"],
		dna_bulge=config["DNAbb"],
		max_mismatch=config["max_mismatch"],
		casoff_accelerator=config["PU"]
	conda:
		"../envs/casoff.yaml"
	message:
		"""
		"""
	script:
		"py/build_casoff_output.py"




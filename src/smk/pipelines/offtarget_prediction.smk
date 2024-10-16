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
		# Decompress the main reference genome for CasOffinder
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{reference_id}.fa",
			fasta_root_path=config["fasta_root_path"], root_dir=config["output_directory"],
			mode=config["processing_mode"], run_name=config["run_name"],
			offtarget_genomes=config["reference_id"],
			reference_id=config["reference_id"]),
		# Create symlinks of consensus fasta files of alternate genomes for CasOffinder
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{offtarget_genomes}.fa",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"],
			offtarget_genomes=config["offtarget_extended"], reference_id=config["reference_id"]),
		# Prepare input files for casoffinder on a per-editor basis
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/input_files/{query_index}_{editing_tool}_casoff_in.txt",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"], reference_id=config["reference_id"],
			offtarget_genomes=config["offtarget_genomes"],
			editing_tool=config["editors_list"],
			query_index=config['query_index']),
		# Run Cas-Offinder
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_casoff.txt",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"],reference_id=config["reference_id"],
			offtarget_genomes=config["offtarget_genomes"],
			editing_tool=config["editors_list"],
			query_index=config['query_index']),
		expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_casoff_parsed.txt",
			root_dir=config["output_directory"],mode=config["processing_mode"],
			run_name=config["run_name"],reference_id=config["reference_id"],
			offtarget_genomes=config["offtarget_genomes"],
			editing_tool=config["editors_list"],
			query_index=config['query_index']),

# noinspection SmkAvoidTabWhitespace
rule decompress_genome:
	input:
		assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{{reference_id}}.fa.gz".format(
			fasta_root_path=config["fasta_root_path"]))
	output:
		decompressed_assembly_symlink = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{reference_id}.fa"
	params:
		decompressed_assembly_path = lambda wildcards: glob.glob("{fasta_root_path}/{reference_id}.fa".format(
			fasta_root_path=config["fasta_root_path"],reference_id=wildcards.reference_id)),
		link_directory = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/"
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
rule symlink_genomes:
	input:
		consensus_fasta = lambda wildcards: glob.glob("{meditdb_path}/{{mode}}/consensus_refs/{{reference_id}}/{{offtarget_genomes}}.fa".format(
			meditdb_path=config["meditdb_path"]))
	output:
		decompressed_assembly_symlink = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{offtarget_genomes}.fa",
	params:
		link_directory = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/",
	shell:
		"""
		ln --symbolic -t {params.link_directory} {input.consensus_fasta}
		"""

# noinspection SmkAvoidTabWhitespace
rule casoff_input_formatting:
	input:
		guides_per_editor_path = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}.pkl",
		guide_search_params = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/dynamic_params/{query_index}_guide_search_params.pkl",
		decompressed_assembly_symlink = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{offtarget_genomes}.fa",
	output:
		casoff_input = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/input_files/{query_index}_{editing_tool}_casoff_in.txt",
		casoff_support = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/input_files/{query_index}_{editing_tool}_casoff_support.txt"
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
		casoff_input="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/input_files/{query_index}_{editing_tool}_casoff_in.txt"
	output:
		casoff_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_casoff.txt"
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
rule casoff_scoring:
	input:
		casoff_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_casoff.txt",
		casoff_support = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/input_files/{query_index}_{editing_tool}_casoff_support.txt",
	output:
		offtarget_scores = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_casoff_scores.pkl",
		formatted_casoff_temp = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_casoff_temp_parsed.txt"
	params:
		rna_bulge=config["RNAbb"],
		dna_bulge=config["DNAbb"],
		max_mismatch=config["max_mismatch"],
		casoff_accelerator=config["PU"],
		annote_path = config["refseq_table"],
		models_path = config["models_path"]
	conda:
		"../envs/casoff.yaml"
	message:
		"""
# === PROCESS OFFTARGET SCORING === #
Inputs used:
--> Analyze off-target effect for guides predicted for: {wildcards.editing_tool}
--> Take formatted inputs from :\n {input.casoff_out}

Run parameters:
--> RNA bulge: {params.rna_bulge} 
--> DNA bulge: {params.dna_bulge}
--> Maximum mismatch: {params.max_mismatch}
--> Cas-Offinder running on device: {params.casoff_accelerator}
--> RefSeq Table: {params.annote_path}
--> Path to pickled models: {params.models_path}

Outputs generated:
--> Scores: {output.offtarget_scores}
--> Intermediate file: {output.formatted_casoff_temp}
Wildcards in this rule:
--> {wildcards}		
		"""
	script:
		"py/build_casoff_scores.py"

# noinspection SmkAvoidTabWhitespace
rule casoff_output_formatting:
	input:
		offtarget_scores = expand("{{root_dir}}/{{mode}}/jobs/{{run_name}}/guide_prediction-{{offtarget_genomes}}/offtarget_prediction/{offtarget_genomes}/{{query_index}}_{editing_tool}_casoff_scores.pkl",
			editing_tool=config['editors_list']),
		formatted_casoff_temp = expand("{{root_dir}}/{{mode}}/jobs/{{run_name}}/guide_prediction-{{offtarget_genomes}}/offtarget_prediction/{offtarget_genomes}/{{query_index}}_{editing_tool}_casoff_temp_parsed.txt",
			editing_tool=config['editors_list'])
	output:
		formatted_casoff_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_Guides_found_casoff.txt",
		offtarget_scores_out = "{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_casoff_scores.txt",
	params:
		rna_bulge=config["RNAbb"],
		dna_bulge=config["DNAbb"],
		max_mismatch=config["max_mismatch"],
	conda:
		"../envs/casoff.yaml"
	message:
		"""
# === COMPILE/FORMAT OFFTARGET OUTPUTS === #
Inputs used:
--> Take formatted inputs from :\n {input.offtarget_scores}

Run parameters:
--> RNA bulge: {params.rna_bulge} 
--> DNA bulge: {params.dna_bulge}
--> Maximum mismatch: {params.max_mismatch}


Outputs generated:
--> Updated Guides_found with CasOff scores: {output.formatted_casoff_out}
--> CasOff scores: {output.offtarget_scores_out}
Wildcards in this rule:
--> {wildcards}				
		"""
	script:
		"py/build_casoff_scores.py"

# noinspection SmkAvoidTabWhitespace
# TODO: Compile all Offtargets reports
rule aggregate_altgenome_reports:
	input:
		a = ""
	output:
		b = ""
	params:
		c = ""
	conda:
		"../envs/vcf.yaml"
	message:
		"""
		"""
	script:
		"py/"

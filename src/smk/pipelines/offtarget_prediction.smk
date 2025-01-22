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
        # Create symlinks of consensus fasta files of alternate genomes for CasOffinder
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{offtarget_genomes}.fa.index",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            run_name=config["run_name"],
            offtarget_genomes=config["offtarget_extended"],reference_id=config["reference_id"]),
        # Prepare input files for casoffinder on a per-editor basis
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/input_files/{query_index}_{editing_tool}_guidescan_in.csv",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            run_name=config["run_name"],reference_id=config["reference_id"],
            offtarget_genomes=config["offtarget_genomes"],
            genome_type=config["genome_type"],
            editing_tool=config["editors_list"],
            query_index=config['query_index'],
            editor_pam=config['pams_list']),
        # Run GuideScan2
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_guidescan_filtered.bed",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            run_name=config["run_name"],reference_id=config["reference_id"],
            offtarget_genomes=config["offtarget_genomes"],
            offtarget_bed_files=config["offtarget_bed_files"],
            editing_tool=config["editors_list"],
            query_index=config['query_index'],
            alt_pam_list=config['alt_pams_list'],
            pam_is_first=config['pam_is_first_list']),
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_guidescan.csv",
            root_dir=config["output_directory"],
            mode=config["processing_mode"],
            run_name=config["run_name"],
            reference_id=config["reference_id"],
            offtarget_genomes=config["offtarget_genomes"],
            editing_tool=config["editors_list"],
            query_index=config['query_index']),


# noinspection SmkAvoidTabWhitespace
rule symlink_genomes:
    input:
        consensus_fasta=lambda wildcards: glob.glob("{meditdb_path}/{{mode}}/consensus_refs/{{reference_id}}/{{offtarget_genomes}}.fa".format(
            meditdb_path=config["meditdb_path"]))
    output:
        decompressed_assembly_symlink="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{offtarget_genomes}.fa.index",
    params:
        link_directory="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/",
    shell:
        """
		ln --symbolic -t {params.link_directory} {input.consensus_fasta}
		"""

# noinspection SmkAvoidTabWhitespace
rule casoff_input_formatting:
    input:
        guides_per_editor_path="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/dynamic_params/{query_index}_{editing_tool}.pkl",
    output:
        casoff_input="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/input_files/{query_index}_{editing_tool}_guidescan_input.csv",
    conda:
        "../envs/casoff.yaml"
    message:
        """
# === DATA FORMATTING FOR Guidescan2 === #	
Inputs used:
--> Take guides grouped by editing tool:\n {input.guides_per_editor_path}
--> Use indexed fasta:\n {input.decompressed_assembly_symlink} 

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
        casoff_input="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/input_files/{query_index}_{editing_tool}_guidescan_input.csv",
        ref_guidescan_full_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{reference_id}/{query_index}_{editing_tool}_guidescan_filtered.bed",#only needed for aextended genomes. If easier make empty bed file for ref
        vcf_bed="/standard/bed_files/{offtarget_genomes}.bed" ## Only needed for an extended genomes. If easier make empty bed file for ref
    output:
        guidescan_tmp_full_csv="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_guidescan_tmp_full.csv",
        guidescan_tmp_full_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_guidescan.bed",
        guidescan_tmp_missing_from_ref_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_missing_from_ref.bed",
        guidescan_tmp_variants_combined_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_combined.bed",
        guidescan_filtered_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_guidescan_filtered.bed"
    params:
        rna_bulge=config["RNAbb"],
        dna_bulge=config["DNAbb"],
        max_mismatch=config["max_mismatch"]
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

Outputs generated:
--> CasOffinder output: {output.casoff_out}
Wildcards in this rule:
--> {wildcards}		
		"""
    shell:
        """
        #### STEP 1 - run guidescan
        
        #PAM is 3' and there are no alt pams. ex: spCas9 NGG
		if [ {wildcard.alt_pam_list} == "no_alt_pam" ] && [ {wildcards.pam_is_first} != "--start"]; then 
		    guidescan enumerate -m {params.max_mismatches} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} -f {input.casoff_input} -n {params.threads} -o {output.guidescan_tmp_full_csv} {fasta}'
		    
		#PAM is 5' and there are alt pams. ex: Cas12a TTTV
		elif [ {wildcard.alt_pam_list} != "no_alt_pam" ] && [ {wildcards.pam_is_first} == "--start"]; then
		    guidescan enumerate -m {params.max_mismatches} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} {wildcards.pam_is_first} --alt-pam {wildcards.alt_pam_list} -f {input.casoff_input} -n {params.threads} -o {output.guidescan_tmp_full_csv} {fasta}'
        
        #PAM is 5' and there are no alt pams. ex: CasX TTCN
        elif [ {wildcard.alt_pam_list} == "no_alt_pam" ]; then
            guidescan enumerate -m {params.max_mismatches} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} {wildcards.pam_is_first} -f {input.casoff_input} -n {params.threads} -o {output.guidescan_tmp_full_csv} {fasta}'
        #PAM is 5' and there are  alt pams. ex: saCas9 NNGRR
        else
            guidescan enumerate -m {params.max_mismatches} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} --alt-pam {wildcards.alt_pam_list} -f {input.casoff_input} -n {params.threads} -o {output.guidescan_tmp_full_csv} {fasta}'
        fi
        
        #### STEP 2: If extended genome. Only keep the sites that differ from ref genome
        
        #Extended 
        if [[ {input.genome_type} == "extended" ]]; then
        
            #convert guidescan csv to bed
            awk -F',' 'NR>1 {print $3 "\t" $4 "\t" $4+30 "\t" $0",added"}' {output.guidescan_tmp_full_csv} | bedtools sort -i > {output.guidescan_tmp_full_bed}
            
            #subset ref sites that are missing the alt guidescan
            bedtools subtract -a {input.ref_guidescan_full_bed} -b {output.guidescan_tmp_full_bed} -wa > {output.guidescan_tmp_missing_from_ref_bed}
    
            #combine the missing ref sites with the full output
            cat {output.guidescan_tmp_full_bed} {output.guidescan_tmp_missing_from_ref_bed} | bedtools sort -i > {output.guidescan_tmp_variants_combined_bed}
            
            #drop anything wihout an overlapping variant (keep only sites though effected by alt variants we don't want the same sites as the reference
            bedtools intersect -a {output.guidescan_tmp_variants_combined_bed} -b {input.vcf_bed} -wa -wb > {output.guidescan_filtered_bed}
            
            rm {output.guidescan_tmp_full_bed}
            rm {output.guidescan_tmp_missing_from_ref_bed}
            rm {output.guidescan_tmp_variants_combined_bed}
        
        
        # Reference -- simply make bed, no filtering needed
        else
            awk -F',' 'NR>1 {print $3 "\t" $4 "\t" $4+30 "\t" $0",removed"}' {output.guidescan_tmp_full_csv} | bedtools sort -i > {output.guidescan_filtered_bed}
        fi
    
        rm {output.guidescan_tmp_full_csv}

       
		"""

# noinspection SmkAvoidTabWhitespace
# ref_out and genome_type needed to compare alt off_targets to reference output
rule casoff_scoring:
    input:
        #Temp file
        guidescan_filtered_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_guidescan_filtered.bed",
    output:
        formatted_casoff="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_Offtargets_found.csv"
    params:
        #genome_type = reference or extended
        rna_bulge=config["RNAbb"],
        dna_bulge=config["DNAbb"],
        max_mismatch=config["max_mismatch"],
    conda:
        "../envs/casoff.yaml"
    message:
        """
# === PROCESS OFFTARGET SCORING === #
Inputs used:
--> Analyze off-target effect for guides predicted for: {wildcards.editing_tool}
--> Take formatted inputs from :\n {input.guidescan_filtered_bed}

Run parameters:
--> RNA bulge: {params.rna_bulge} 
--> DNA bulge: {params.dna_bulge}
--> RefSeq Table: {params.annote_path}
--> Path to pickled models: {params.models_path}

Outputs generated:
--> Reformatted Guide scan file: {output.formatted_casoff}
Wildcards in this rule:
--> {wildcards}		
		"""
    script:
        "py/build_casoff_scores.py"
#TODO: remove guidescan_filtered_bed, it's a temperay file used to make formatted_casoff

# === COMPILE/FORMAT OFFTARGET OUTPUTS === #
# noinspection SmkAvoidTabWhitespace
rule casoff_output_formatting:
    input:
        off_target_directory = ""
        # just need directory
    output:
        offtarget_summary="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/summary_reports/{query_index}_offtarget_summary.csv",
        off_target_summary_expanded="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/summary_reports/{query_index}_offtarget_summary_expanded.csv"
    params:
        all_offtarget_genomes=config["offtarget_genomes"],
        editors_list=config["editor_list"],
        rna_bulge=config["RNAbb"],
        dna_bulge=config["DNAbb"],
        max_mismatch=config["max_mismatch"],
    conda:
        "../envs/casoff.yaml"
    message:
        """
# === COMPILE/FORMAT OFFTARGET OUTPUTS === #
Inputs used:
--> Take formatted inputs from :\n {input.off_target_directory}

Run parameters:
--> RNA bulge: {params.rna_bulge} 
--> DNA bulge: {params.dna_bulge}
--> Maximum mismatch: {params.max_mismatch}
--> list of editors: {params.editors_list}

Outputs generated:
--> Aggregate summary of all genome off-target sites: {output.offtarget_summary}
--> Expanded version of aggregate all genome off-target sites: {output.offtarget_summary_expanded}
Wildcards in this rule:
--> {wildcards}				
		"""
    script:
        "py/build_casoff_output.py"


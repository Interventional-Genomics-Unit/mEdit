# **** Variables ****
configfile: ""

import glob

# noinspection SmkAvoidTabWhitespace
rule all:
    input:
        # === Create consensus FASTA sequences for Alternative Genomes
        #   == rule
        expand("{meditdb_path}/{mode}/consensus_refs/{reference_id}/{offtarget_extended}.fa",
            meditdb_path=config["meditdb_path"],mode=config["processing_mode"],
            reference_id=config["sequence_id"],offtarget_extended=config["offtarget_extended"]),
        # === Setup Guide Scan Indices ===
        #   == rule set_gscan_indices ==
        expand("{root_dir}/{mode}/tmp/{offtarget_extended}.consensus.fa.index.gs",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            offtarget_extended=config["offtarget_extended"]),
        # === Prepare input files for casoffinder on a per-editor basis ===
        #   == rule casoff_input_formatting ==
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/input_files/{query_index}_{editing_tool}_guidescan_input.csv",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            run_name=config["run_name"],reference_id=config["sequence_id"],
            query_index=config['query_index'],editing_tool=config["pam_per_editor_dict"],
            ),
        # === Run GuideScan2 ===
        #   == rule casoff_run_ref ==
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{reference_id}/{query_index}_{editing_tool}_guidescan_filtered.bed",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            run_name=config["run_name"],reference_id=config["sequence_id"],
            query_index=config['query_index'],editing_tool=config["pam_per_editor_dict"]),
        #   == rule casoff_run_extended ==
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_extended}/{query_index}_{editing_tool}_guidescan_filtered.bed",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            run_name=config["run_name"],reference_id=config["sequence_id"],
            offtarget_extended=config["offtarget_extended"],editing_tool=config["pam_per_editor_dict"],
            query_index=config['query_index']),

        # === Post-rprocessing and Formatting of Guidescam Output Files ===
        #   == rule casoff_scoring ==
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_Offtargets_found.csv",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            run_name=config["run_name"],reference_id=config["sequence_id"],
            offtarget_genomes=config["offtarget_genomes"],query_index=config['query_index'],
            editing_tool=config["pam_per_editor_dict"]),
        #   == rule casoff_output_formatting ==
        expand("{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/summary_reports/{query_index}_offtarget_summary.csv",
            root_dir=config["output_directory"],mode=config["processing_mode"],
            run_name=config["run_name"],reference_id=config["sequence_id"],
            query_index=config['query_index'])

# noinspection SmkAvoidTabWhitespace
rule set_consensus_fasta:
    input:
        assembly_path=lambda wildcards: glob.glob("{fasta_root_path}/{reference_id}.filtered.fa.gz".format(
            fasta_root_path=config["fasta_root_path"],reference_id=wildcards.reference_id)),
        filtered_vcf="{meditdb_path}/{mode}/consensus_refs/{reference_id}/{offtarget_extended}.filtered.vcf.gz"
    output:
        consensus_fasta="{meditdb_path}/{mode}/consensus_refs/{reference_id}/{offtarget_extended}.fa",
        bed_file="{meditdb_path}/{mode}/bed_files/{reference_id}/{offtarget_extended}.bed"
    params:
        tmp_vcf = "{meditdb_path}/{mode}/consensus_refs"
    conda:
        "../envs/samtools.yaml"
    message:
        """
# === CREATING CONSENSUS FASTA FOR GENOME {wildcards.offtarget_extended} === #	
Input used:
--> Take filtered VCF input file from:\n {input.filtered_vcf}
Output generated:
--> Export Consensus sequence to:\n {output.consensus_fasta}
        """
    shell:
        """
bcftools query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' {input.filtered_vcf} > {output.bed_file}

#bcftools will skip Multiallelic positions when creatign a consensus fasta.
#This removes one of the alt alleles which can be reffered back to in the bedfile later in the pipeline
gunzip -c {input.filtered_vcf} | awk 'BEGIN {{OFS="\t"}} /^#/ {{print; next}} !seen[$1":"$2]++' | bgzip > {params.tmp_vcf}/tmp_{wildcards.offtarget_extended}.vcf
bcftools index -f -t {params.tmp_vcf}/tmp_{wildcards.offtarget_extended}.vcf

## Consensus Sequence
# remove 0|0, .|0, 0| and then chose the alt allele to make fasta
samtools dict {input.assembly_path}
samtools faidx {input.assembly_path}	

bcftools consensus -H A -M A --fasta-ref {input.assembly_path} {params.tmp_vcf}/tmp_{wildcards.offtarget_extended}.vcf > {output.consensus_fasta}
samtools faidx {output.consensus_fasta}        	
        """

# noinspection SmkAvoidTabWhitespace
rule set_gscan_indices:
    input:
        consensus_fasta=lambda wildcards: glob.glob("{meditdb_path}/{mode}/consensus_refs/{reference_id}/{offtarget_extended}.fa".format(
            meditdb_path=config["meditdb_path"], mode=wildcards.mode,
            reference_id=config["sequence_id"], offtarget_extended=wildcards.offtarget_extended))
    output:
        gs_index="{root_dir}/{mode}/tmp/{offtarget_extended}.consensus.fa.index.gs"
    params:
        gscan_indices_dir=lambda wildcards: glob.glob("{meditdb_path}/gscan_indices".format(
            meditdb_path=config["meditdb_path"])),
        gs_index_prefix=lambda wildcards: glob.glob("{meditdb_path}/gscan_indices/{offtarget_genomes}.consensus.fa.index".format(
            meditdb_path=config["meditdb_path"], offtarget_genomes=wildcards.offtarget_extended)),
        gs_index_dummy_dir="{root_dir}/{mode}/tmp"
    conda:
        "../envs/gscan.yaml"
    message:
        """
# === CREATING GuideScan 2 INDEX FOR GENOME {wildcards.offtarget_extended} === #	
Input used:
--> Take Guidescan input file from:\n {input.consensus_fasta}
Output generated:
--> Export Guidescan index to:\n {output.gs_index}
        """
    shell:
        """
cd {params.gscan_indices_dir}        
guidescan index --index {wildcards.offtarget_extended}.consensus.fa.index {input.consensus_fasta}
ln --symbolic -f {params.gscan_indices_dir}/{wildcards.offtarget_extended}.consensus.fa.index.forward {params.gs_index_dummy_dir}/{wildcards.offtarget_extended}.consensus.fa.index.forward
ln --symbolic -f {params.gscan_indices_dir}/{wildcards.offtarget_extended}.consensus.fa.index.reverse {params.gs_index_dummy_dir}/{wildcards.offtarget_extended}.consensus.fa.index.reverse
ln --symbolic -f {params.gscan_indices_dir}/{wildcards.offtarget_extended}.consensus.fa.index.gs {params.gs_index_dummy_dir}/{wildcards.offtarget_extended}.consensus.fa.index.gs
        """

# noinspection SmkAvoidTabWhitespace
rule casoff_input_formatting:
    input:
        guides_per_editor_path="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/dynamic_params/{query_index}_{editing_tool}.pkl"
    output:
        casoff_input="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/input_files/{query_index}_{editing_tool}_guidescan_input.csv"
    params:
        pam_per_editor_dict=lambda wildcards: config["pam_per_editor_dict"][wildcards.editing_tool]
    conda:
        "../envs/gscan.yaml"
    message:
        """
# === INPUT FORMATTING FOR Guidescan 2 [REFERENCE GENOME] === #	
Editing tool processed in this job: {wildcards.editing_tool} 
Inputs used:
--> Take guides grouped by editing tool:\n {input.guides_per_editor_path}
Outputs generated:
--> Guidescan formatted input: {output.casoff_input}
Wildcards in this rule:
--> {wildcards}
"""
    script:
        "py/build_casoff_input.py"

# noinspection SmkAvoidTabWhitespace
rule casoff_run_ref:
    input:
        casoff_input="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/input_files/{query_index}_{editing_tool}_guidescan_input.csv",
        gs_index=lambda wildcards: glob.glob(f"{config['meditdb_path']}/gscan_indices/{wildcards.reference_id}.consensus.fa.index.gs")
    output:
        ref_guidescan_full_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{reference_id}/{query_index}_{editing_tool}_guidescan_filtered.bed"
    params:
        # == Temp Files
        ref_guidescan_full_csv="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{reference_id}/{query_index}_{editing_tool}_guidescan_filtered.csv",
        # == GuideScan Indices Path
        gscan_indices_path=lambda wildcards: glob.glob(f"{config['meditdb_path']}/gscan_indices"),
        # == GuideScan Parameters
        pam_is_first=lambda wildcards: config["pam_is_first_per_editor_dict"][wildcards.reference_id][
            wildcards.editing_tool],
        alt_pams=lambda wildcards: config["alt_pam_per_editor_dict"][wildcards.reference_id][
            wildcards.editing_tool],
        rna_bulge=config["RNAbb"],
        dna_bulge=config["DNAbb"],
        max_mismatch=config["max_mismatch"]
    conda:
        "../envs/gscan.yaml"
    threads:
        int(config["threads"])
    message:
        """
# === PREDICT OFFTARGET EFFECT [REFERENCE GENOME] === #	
Inputs used:
--> Take Guidescan input file from:\n {input.casoff_input}
--> Guide Scan Index:\n {input.gs_index}
--> Filtered Guide Scan csv:\n {params.ref_guidescan_full_csv}
Outputs generated:
--> Generate Reference Genome Bed file: {output.ref_guidescan_full_bed}
Wildcards in this rule:
--> {wildcards}
"""
    shell:
        """
        ## == REFERENCE GENOMES ==
        #### STEP 1 - run guidescan
        
        #PAM is 3' and there are no alt pams. ex: spCas9 NGG
		if [ {params.alt_pams} == "no_alt_pam" ] && [ {params.pam_is_first} != "--start"]; then 
		    guidescan enumerate {params.gscan_indices_path}/{wildcards.reference_id}.consensus.fa.index -m {params.max_mismatch} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} -f {input.casoff_input} -n {threads} -o {params.ref_guidescan_full_csv}
		    
		#PAM is 5' and there are alt pams. ex: Cas12a TTTV
		elif [ {params.alt_pams} != "no_alt_pam" ] && [ {params.pam_is_first} == "--start"]; then
		    guidescan enumerate {params.gscan_indices_path}/{wildcards.reference_id}.consensus.fa.index -m {params.max_mismatch} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} {params.pam_is_first} --alt-pam {params.alt_pams} -f {input.casoff_input} -n {threads} -o {params.ref_guidescan_full_csv}
        
        #PAM is 5' and there are no alt pams. ex: CasX TTCN
        elif [ {params.alt_pams} == "no_alt_pam" ]; then
            guidescan enumerate {params.gscan_indices_path}/{wildcards.reference_id}.consensus.fa.index -m {params.max_mismatch} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} {params.pam_is_first} -f {input.casoff_input} -n {threads} -o {params.ref_guidescan_full_csv}
        #PAM is 5' and there are  alt pams. ex: saCas9 NNGRR
        else
            guidescan enumerate {params.gscan_indices_path}/{wildcards.reference_id}.consensus.fa.index -m {params.max_mismatch} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} --alt-pam {params.alt_pams} -f {input.casoff_input} -n {threads} -o {params.ref_guidescan_full_csv}
        fi
        
        #convert guidescan csv to bed
        awk -F',' 'NR>1 && $4!="NA" {{print $3 "\t" $4-1 "\t" $4+30 "\t" $0",removed"}}' {params.ref_guidescan_full_csv}  | bedtools sort -i > {output.ref_guidescan_full_bed}

        rm {params.ref_guidescan_full_csv}
        """

# noinspection SmkAvoidTabWhitespace
rule casoff_run_extended:
    input:
        casoff_input="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/input_files/{query_index}_{editing_tool}_guidescan_input.csv",
        ref_guidescan_full_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{reference_id}/{query_index}_{editing_tool}_guidescan_filtered.bed",
        bed_file=lambda wildcards: glob.glob("{meditdb_path}/{mode}/bed_files/{reference_id}/{offtarget_extended}.bed".format(
            meditdb_path=config["meditdb_path"], mode=wildcards.mode,reference_id=wildcards.reference_id, offtarget_extended=wildcards.offtarget_extended)),
        gs_index="{root_dir}/{mode}/tmp/{offtarget_extended}.consensus.fa.index.gs",
    output:
        guidescan_filtered_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_extended}/{query_index}_{editing_tool}_guidescan_filtered.bed"
    params:
        # == Temp Files
        guidescan_tmp_full_csv="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_extended}/{query_index}_{editing_tool}_guidescan_tmp_full.csv",
        guidescan_tmp_full_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_extended}/{query_index}_{editing_tool}_guidescan.bed",
        guidescan_tmp_missing_from_ref_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_extended}/{query_index}_{editing_tool}_missing_from_ref.bed",
        guidescan_tmp_variants_combined_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_extended}/{query_index}_{editing_tool}_combined.bed",
        # == GuideScan Indices Path
        gscan_indices_path=lambda wildcards: glob.glob("{meditdb_path}/gscan_indices".format(meditdb_path=config["meditdb_path"])),
        # == GuideScan Parameters
        pam_is_first=lambda wildcards: config["pam_is_first_per_editor_dict"][wildcards.offtarget_extended][wildcards.editing_tool],
        alt_pams=lambda wildcards: config["alt_pam_per_editor_dict"][wildcards.offtarget_extended][wildcards.editing_tool],
        rna_bulge=config["RNAbb"],
        dna_bulge=config["DNAbb"],
        max_mismatch=config["max_mismatch"],
        offtarget_extended=config["offtarget_extended"]
    conda:
        "../envs/gscan.yaml"
    threads:
        int(config["threads"])
    message:
        """
# === PREDICT OFFTARGET EFFECT [ADDITIONAL GENOMES] === #
Inputs used:
--> Analyze off-target effect for guides predicted for: {wildcards.editing_tool}
--> Take formatted inputs from :\n {input.casoff_input}

Run parameters:
--> RNA bulge: {params.rna_bulge} 
--> DNA bulge: {params.dna_bulge}
--> Maximum mismatch: {params.max_mismatch}

Outputs generated:
--> Guidescan output: {output.guidescan_filtered_bed}
Wildcards in this rule:
--> {wildcards}		
		"""
    shell:
        """
        #### STEP 1 - run guidescan        
        #PAM is 3' and there are no alt pams. ex: spCas9 NGG
		if [ {params.alt_pams} == "no_alt_pam" ] && [ {params.pam_is_first} != "--start"]; then 
		    guidescan enumerate {params.gscan_indices_path}/{wildcards.offtarget_extended}.consensus.fa.index -m {params.max_mismatch} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} -f {input.casoff_input} -n {threads} -o {params.guidescan_tmp_full_csv}
		    
		#PAM is 5' and there are alt pams. ex: Cas12a TTTV
		elif [ {params.alt_pams} != "no_alt_pam" ] && [ {params.pam_is_first} == "--start"]; then
		    guidescan enumerate {params.gscan_indices_path}/{wildcards.offtarget_extended}.consensus.fa.index -m {params.max_mismatch} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} {params.pam_is_first} --alt-pam {params.alt_pams} -f {input.casoff_input} -n {threads} -o {params.guidescan_tmp_full_csv}
        
        #PAM is 5' and there are no alt pams. ex: CasX TTCN
        elif [ {params.alt_pams} == "no_alt_pam" ]; then
            guidescan enumerate {params.gscan_indices_path}/{wildcards.offtarget_extended}.consensus.fa.index -m {params.max_mismatch} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} {params.pam_is_first} -f {input.casoff_input} -n {threads} -o {params.guidescan_tmp_full_csv}
        #PAM is 5' and there are  alt pams. ex: saCas9 NNGRR
        else
            guidescan enumerate {params.gscan_indices_path}/{wildcards.offtarget_extended}.consensus.fa.index -m {params.max_mismatch} --rna-bulges {params.rna_bulge} --dna-bulges {params.dna_bulge} --alt-pam {params.alt_pams} -f {input.casoff_input} -n {threads} -o {params.guidescan_tmp_full_csv}
        fi
        
        #### STEP 2: Only keep the sites that differ from ref genome        
        #convert guidescan csv to bed
        awk -F',' 'NR>1 && $4!="NA" {{print $3 "\t" $4-1 "\t" $4+30 "\t" $0",added"}}' {params.guidescan_tmp_full_csv}  | bedtools sort -i > {params.guidescan_tmp_full_bed}
        
        #subset ref sites that are missing the alt guidescan
        bedtools subtract -a {input.ref_guidescan_full_bed} -b {params.guidescan_tmp_full_bed} -wa > {params.guidescan_tmp_missing_from_ref_bed}

        #combine the missing ref sites with the full output
        cat {params.guidescan_tmp_full_bed} {params.guidescan_tmp_missing_from_ref_bed} | bedtools sort -i > {params.guidescan_tmp_variants_combined_bed}
        
        #drop anything wihout an overlapping variant (keep only sites though effected by alt variants we don't want the same sites as the reference
        bedtools intersect -a {params.guidescan_tmp_variants_combined_bed} -b {input.bed_file} -wa -wb > {output.guidescan_filtered_bed}
        
        rm {params.guidescan_tmp_full_bed}
        rm {params.guidescan_tmp_missing_from_ref_bed}
        rm {params.guidescan_tmp_variants_combined_bed}    
        rm {params.guidescan_tmp_full_csv}
		"""

# noinspection SmkAvoidTabWhitespace
rule casoff_scoring:
    input:
        #Temp file
        guidescan_filtered_bed="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_guidescan_filtered.bed"
    output:
        formatted_casoff="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/{offtarget_genomes}/{query_index}_{editing_tool}_Offtargets_found.csv"
    params:
        extended_genomes=config["offtarget_extended"],
        rna_bulge=config["RNAbb"],
        dna_bulge=config["DNAbb"],
        models_path=config["models_path"],
        annote_path=config["refseq_table"],
    conda:
        "../envs/medit.yaml"
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

# noinspection SmkAvoidTabWhitespace
rule casoff_output_formatting:
    input:
        guides_per_editor_list=expand("{{root_dir}}/{{mode}}/jobs/{{run_name}}/guide_prediction-{{reference_id}}/offtarget_prediction/dynamic_params/{{query_index}}_{editing_tool}.pkl",
            editing_tool=config["pam_per_editor_dict"]),
        formatted_casoff_list=expand("{{root_dir}}/{{mode}}/jobs/{{run_name}}/guide_prediction-{{reference_id}}/offtarget_prediction/{offtarget_genomes}/{{query_index}}_{editing_tool}_Offtargets_found.csv",
            offtarget_genomes=config["offtarget_genomes"],editing_tool=config["pam_per_editor_dict"])
    output:
        offtarget_summary="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/summary_reports/{query_index}_offtarget_summary.csv",
        off_target_summary_expanded="{root_dir}/{mode}/jobs/{run_name}/guide_prediction-{reference_id}/offtarget_prediction/summary_reports/{query_index}_offtarget_summary_expanded.csv"
    params:
        extended_genomes=config["offtarget_extended"],
        max_mismatch=config["max_mismatch"],
        rna_bulge=config["RNAbb"],
        dna_bulge=config["DNAbb"]
    conda:
        "../envs/gscan.yaml"
    message:
        """
# === COMPILE/FORMAT OFFTARGET OUTPUTS === #
Inputs used:
--> Take formatted inputs from :\n {input.formatted_casoff_list}

Run parameters:
--> RNA bulge: {params.rna_bulge} 
--> DNA bulge: {params.dna_bulge}
--> Maximum mismatch: {params.max_mismatch}

Outputs generated:
--> Aggregate summary of all genome off-target sites: {output.offtarget_summary}
--> Expanded version of aggregate all genome off-target sites: {output.off_target_summary_expanded}
Wildcards in this rule:
--> {wildcards}				
		"""
    script:
        "py/build_casoff_output.py"

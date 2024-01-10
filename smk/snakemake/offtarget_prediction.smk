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
rule:
	input:
		guides_report_out = ""
	output:
		casoff_out = ""
	params:
		tmp_casoff = ""
	conda:
	threads:
	message:
	shell:

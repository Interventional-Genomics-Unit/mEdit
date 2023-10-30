#!/bin/bash

## Fetching, cleaning and creation of a Fasta file from a VCF file
####

OUTDIR=$1
PREFIX=$2 #HG02257
REF=$3 #hg38_standard.fa.gz


# 1) download diploid VCF files from AWS (-->to be a loop using index file) 
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/"${PREFIX}"/assemblies/year1_freeze_assembly_v2/assembly_qc/dipcall_v0.2/"${PREFIX}".f1_assembly_v2.dip.vcf.gz

# 2) filter GAT1 or GAT2 samples (samples where one haplotype has a sequence depth = 0)
# & filter reference & variant alleles > 5nt
# Create index file
bcftools filter -O z -o "${PREFIX}".filtered.vcf.gz -e 'GT="." || ILEN <= -5 || ILEN >= 5' "${VCF}" | bcftools index -t "${PREFIX}".f1_assembly_v2.dip.vcf.gz -i

# 3) Making a consensus
#previously made a seperate hg38 Ref Fasta that only have standard chromsomes --> /groups/clinical/projects/editability/tables/raw_tables/VCFs/hg38_standard.fa.gz

samtools dict "${REF}" -o "${REF}".dict
# shellcheck disable=SC2086
gatk FastaAlternateReferenceMaker -R ${REF} -O "${PREFIX}".consensus.fa.gz -V "${PREFIX}".filtered.vcf.gz --sequence-dictionary "${REF}".dict

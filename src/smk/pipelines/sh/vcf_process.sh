#!/usr/bin/env bash

set -euo pipefail

# Assign inputs from Snakemake to variables
source_vcf="${snakemake_input[source_vcf]}"
assembly_path="${snakemake_input[assembly_path]}"
meditdb_path="${snakemake_wildcards[meditdb_path]}"
mode="${snakemake_wildcards[mode]}"
vcf_id="${snakemake_wildcards[vcf_id]}"
filtered_vcf="${snakemake_output[filtered_vcf]}"

# Step 0: Ensure output directory exists
mkdir -p "$(dirname "$filtered_vcf")"

if [ "$mode" == 'vcf' ]; then
  # Derived file paths
  echo "processing vcf mode - vcf"
  decompressed_vcf="${meditdb_path}/${mode}/source_vcfs/${vcf_id}.vcf"
  bgzipped_vcf="${meditdb_path}/${mode}/source_vcfs/bgz_${vcf_id}.vcf.gz"

  mkdir -p "$(dirname "$bgzipped_vcf")"

  # Step 1: Decompress VCF (input files are expected either in VCF, or gzip compressed VCF)
  #   -> This is handled at the 'prog' level, which ensures only Gzip-compressed files are stored by the user
  echo "decompressing '$source_vcf'"
  gunzip -kf "$source_vcf"

  # Step 2: Compress with bgzip
  #   -> Bgzip is a requirement for bcftools;
  #   -> This is handled internally to remove one potential entry barrier for used: handling the Tabix Dependency
  echo "making bgzip file '$bgzipped_vcf'"
  bgzip "$decompressed_vcf" -o "$bgzipped_vcf"
  bcftools index "$bgzipped_vcf"

else
  # Step 2.5: Standard Mode can bypass the Bgzip routine
  #   -> The built-in VCF is already Bgzipped
  echo "processing vcf mode - standard hprc"
  bgzipped_vcf="$source_vcf"
fi

# Step 5: Get sample and contig info
num_samples=$(bcftools query -l "$bgzipped_vcf" | wc -l)
CONTIGS=$(cut -f1 "${assembly_path}.fai" | paste -sd, -)
echo "Number of VCF samples = '$num_samples'"

# Step 6: Apply conditional filters
if [ "$num_samples" -le 1 ]; then
    if bcftools view -h "$bgzipped_vcf" | grep -q '##FORMAT=<ID=DP'; then
        echo "filtering QUAL<15 || FMT/DP<5, extracting only SNVs and splitting multi-allelic per line"
        bcftools view -e 'QUAL<15 || FMT/DP<5' -r "$CONTIGS" "$bgzipped_vcf" | \
        bcftools norm --multiallelics -any -f "$assembly_path" | \
        bcftools view --min-ac=1 -v snps -e 'GT="."' | \
        bcftools sort -W=tbi -O z -o "$filtered_vcf"
    else
        echo "extracting only SNVs and splitting multi-allelic per line"
        bcftools view -r "$CONTIGS" "$bgzipped_vcf" | \
        bcftools norm --multiallelics -any -f "$assembly_path" | \
        bcftools view --min-ac=1 -v snps -e 'GT="."' | \
        bcftools sort -W=tbi -O z -o "$filtered_vcf"
    fi
  if [ "$mode" == 'vcf' ]; then
    # Step 7: Cleanup temporary files on vcf mode
    echo "removing  '$decompressed_vcf' '$bgzipped_vcf' '$bgzipped_vcf.*'"
    rm -f "$decompressed_vcf" "$bgzipped_vcf" "$bgzipped_vcf.*"
  fi
else
  if [ "$mode" == 'vcf' ]; then
    echo "The VCF file '$source_vcf' contains more than one genome. Please provide a single-sample VCF." >&2
    exit 1

  else
    # Standard Mode will land here.
    bcftools view -s "$vcf_id" --min-ac=1 "$bgzipped_vcf" | bcftools sort -W=tbi -O z -o "$filtered_vcf"
  fi
fi

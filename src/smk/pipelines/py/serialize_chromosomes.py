# == Native Modules
import pickle
import re
import subprocess
import sys
# == Installed Modules
from Bio import SeqIO
# == Project Modules


def pickle_chromosomes(genome_fasta, output_dir):
	manifest = []
	records = SeqIO.parse(open(genome_fasta, 'rt'), "fasta")
	for record in records:
		if re.search(r"chr\w{0,2}$", record.id, re.IGNORECASE):
			manifest.append(record.id)
			outfile = f"{output_dir}/{record.id}.pkl"
			with open(outfile, 'ab') as gfile:
				pickle.dump(record, gfile)
	return manifest


def main():
	# === Inputs ===
	assembly_path = str(snakemake.input.assembly_path)
	# === Outputs ===
	serialized_chr_manifest = str(snakemake.output.serialized_chr_manifest)
	# === Params ===
	output_dir = str(snakemake.params.output_dir)
	decompressed_assembly = str(snakemake.params.decompressed_assembly)

	threads = str(snakemake.threads)

	# Cleanup possible whitespaces in the assembly path
	clean_assembly_path = re.sub(r'\s+', '', assembly_path)

	subprocess.run(["bgzip", "-df", "-@", str(threads), clean_assembly_path], check=True, stdout=open(decompressed_assembly, "wb"))

	chr_manifest = pickle_chromosomes(decompressed_assembly, output_dir)

	with open(serialized_chr_manifest, 'w') as file_handle:
		for chromosome in chr_manifest:
			file_handle.write(chromosome)

	subprocess.run(["bgzip", "-cf", "-@", str(threads), decompressed_assembly], check=True,
				   stdout=open(clean_assembly_path, "wb"))


if __name__ == "__main__":
	main()

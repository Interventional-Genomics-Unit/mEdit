# == Native Modules
import pickle
import re
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
	assembly_path = sys.argv[1]  # str(snakemake.input.assembly_path)
	# === Outputs ===
	serialized_chr_manifest = sys.argv[2]  # str(snakemake.output.serialized_chr_manifest)
	# === Params ===
	output_dir = sys.argv[3]  # str(snakemake.params.output_dir)

	chr_manifest = pickle_chromosomes(assembly_path, output_dir)

	with open(serialized_chr_manifest, 'w') as file_handle:
		for chromosome in chr_manifest:
			file_handle.write(chromosome)


if __name__ == "__main__":
	main()

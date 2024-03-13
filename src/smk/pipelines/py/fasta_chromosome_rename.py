# == Native Modules
import re
import sys
import sys
# == Installed Modules
from Bio import SeqIO
# == Project Modules


def main():
	INPUT_FASTA = sys.argv[1]

	record = SeqIO.parse(open(INPUT_FASTA), "fasta")
	changed_recs = []
	for rec in record:
		a = re.search(r"Homo sapiens chromosome (\d+), \S+ Primary Assembly", rec.description)
		try:
			rec.id = f"chr{a.group(1)}"
			rec.name = f"chr{a.group(1)}"
			rec.description = ""
		except AttributeError:
			pass
		changed_recs.append(rec)

		SeqIO.write(changed_recs, INPUT_FASTA, 'fasta')


if __name__ == "__main__":
	main()

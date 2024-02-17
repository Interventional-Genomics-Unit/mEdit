# == Native Modules
import pickle
# == Installed Modules
import pandas as pd
# == Project Modules


def cas_offinder_bulge(input_filename, output_filename,cas_off_expath,bulge):
    '''
     The cas-offinder off-line package contains a bug that doesn't allow bulges
    This script is partially a wrapper for cas-offinder to fix this bug
     created by...
    https://github.com/hyugel/cas-offinder-bulge
    '''
    # INPUT LEG
    fnhead = input_filename.replace("_input.txt", "")
    id_dict = {}
    if bulge == True:
        with open(input_filename) as f:
            chrom_path = f.readline()
            pattern, bulge_dna, bulge_rna = f.readline().strip().split()
            isreversed = False
            for i in range(int(len(pattern) / 2)):
                if pattern[i] == 'N' and pattern[len(pattern) - i - 1] != 'N':
                    isreversed = False
                    break
                elif pattern[i] != 'N' and pattern[len(pattern) - i - 1] == 'N':
                    isreversed = True
                    break
            bulge_dna, bulge_rna = int(bulge_dna), int(bulge_rna)
            targets = [line.strip().split() for line in f]
            rnabulge_dic = defaultdict(lambda: [])
            bg_tgts = defaultdict(lambda: set())
            for raw_target, mismatch, gid in targets:
                if isreversed:
                    target = raw_target.lstrip('N')
                    len_pam = len(raw_target) - len(target)
                    bg_tgts['N' * len_pam + target + 'N' * bulge_dna].add(mismatch)
                    id_dict['N' * len_pam + target + 'N' * bulge_dna] = gid
                    for bulge_size in range(1, bulge_dna+1):
                        for i in range(1, len(target)):
                            bg_tgt = 'N' * len_pam + target[:i] + 'N' * bulge_size + target[i:] + 'N' * (bulge_dna - bulge_size)
                            bg_tgts[bg_tgt].add(mismatch)
                            id_dict[bg_tgt] = gid

                    for bulge_size in range(1, bulge_rna+1):
                        for i in range(1, len(target)-bulge_size):
                            bg_tgt = 'N' * len_pam + target[:i] + target[i+bulge_size:] + 'N' * (bulge_dna + bulge_size)
                            bg_tgts[bg_tgt].add(mismatch)
                            rnabulge_dic[bg_tgt].append((i, int(mismatch), target[i:i+bulge_size]))
                            id_dict[bg_tgt] = gid
                else:
                    target = raw_target.rstrip('N')
                    len_pam = len(raw_target) - len(target)
                    bg_tgts['N' * bulge_dna + target + 'N' * len_pam].add(mismatch)
                    id_dict['N' * bulge_dna + target + 'N' * len_pam] = gid
                    for bulge_size in range(1, bulge_dna+1):
                        for i in range(1, len(target)):
                            bg_tgt = 'N' * (bulge_dna - bulge_size) + target[:i] + 'N' * bulge_size + target[i:] + 'N' * len_pam
                            bg_tgts[bg_tgt].add(mismatch)
                            id_dict[bg_tgt] = gid

                    for bulge_size in range(1, bulge_rna+1):
                        for i in range(1, len(target)-bulge_size):
                            bg_tgt = 'N' * (bulge_dna + bulge_size) + target[:i] + target[i+bulge_size:] + 'N' * len_pam
                            bg_tgts[bg_tgt].add(mismatch)
                            rnabulge_dic[bg_tgt].append( (i, int(mismatch), target[i:i+bulge_size]) )
                            id_dict[bg_tgt] = gid
            if isreversed:
                seq_pam = pattern[:len_pam]
            else:
                seq_pam = pattern[-len_pam:]
        with open(fnhead + '_bulge.txt', 'w') as f:
            f.write(chrom_path)
            if isreversed:
                f.write(pattern + bulge_dna*'N' + '\n')
            else:
                f.write(bulge_dna*'N' + pattern + '\n')
            cnt = 0
            for tgt, mismatch in bg_tgts.items():
                f.write(tgt + ' ' + str(max(mismatch)) + ' ' + '\n')
                cnt+=1
        # THIS FILE PATH IS SUPPLIED TO CASOFF-FINDER
        casin = fnhead + '_bulge.txt'
    else:
        nobulge_dict = {}
        with open(input_filename) as inf:
            for line in inf:
                entry = line.strip().split(' ')
                if len(entry) > 2 and len(entry[-1]) > 3:
                    seq, mm, gid = entry
                    nobulge_dict[seq] = [gid, mm]
        casin = input_filename

    print("Created temporary file (%s)." % (casin))
    # THIS FILE PATH IS SUPPLIED TO CASOFF-FINDER
    outfn = fnhead + '_temp.txt'
    print("Running Cas-OFFinder (output file: %s)..." % outfn)


def run_casoffinder(resultsfolder,
					fasta_fname,
					guide_tab_fname,
					search_params,
					cas_off_expath,
					genome_name,
					guide_src_name,
					casoff_params,
					annote_path):
	#guide_tab_fname = '/groups/clinical/projects/editability/medit_queries/medit_test/test_out/hg38_Guides_found.csv'
	gdf = pd.read_csv(guide_tab_fname)
	ots = {}
	gpr = gdf.groupby('Editor')
	if casoff_params[1:3] == (0, 0):
		bulge = False
	else:
		bulge = True
	# for each editor type find off_targets
	# THIS LOOP AND SUBSEQUENT PANDAS READ/GROUPING/FORMATTING WILL BE INSTANTIATED AT THE prog LEVEL
	for editor, stats in gpr:
		infile = f"{resultsfolder}{genome_name}_{guide_src_name}_{editor}_casoffinder_input.txt"
		pam, pamISfirst, guidelen = search_params[editor][0:3]
		guides, gnames = list(stats.gRNA), list(stats.Guide_ID)
		# make input file
		make_casoffinder_input(infile,
							   fasta_fname,
							   pam,
							   pamISfirst,
							   guidelen,
							   guides,
							   gnames,
							   casoff_params)

		output_filename = infile.replace('_input.txt', '_output.txt')
		# cas_offinder_bulge FUNCTION STARTS HERE
		cas_offinder_bulge(infile, output_filename, cas_off_expath, bulge)


def main():
	'''
	  ### For Daniel to snakemake <---------------

	  input paths for this script:
		  -resultsfolder -- results output folder
		  -guide_search_params -- search paramters used in fetchguides
		  -guide_tab_fname -- original guide table output from FetchGuides OR ALT process_genomes files
		  -fasta_fname -- Hg38 fasta or if using alternative consensus genome
		  -(maybe?) casoffinder path

	  input variables for the script:
		  -genome_name -- name of fasta/consensus we are searching
		  -guides_src_name -- name of the guides source genome ex. HG38 or HG02257

	  *possibly allow for changes in the cas-offinder parameters
	  see bottom of page
	  '''
	# SNAKEMAKE IMPORTS
	# === Inputs ===
	guides_report = str(snakemake.input.guides_report_out)
	fasta_ref = str(snakemake.input.fasta_ref)
	guide_search_params = str(snakemake.input.guide_search_params)
	# snv_site_info = str(snakemake.input.snv_site_info)
	# annote_path = str(snakemake.params.annote_path)
	# === Outputs ===
	casoff_out = str(snakemake.output.casoff_out)
	# === Params ===
	RNAbb = str(snakemake.params.rna_bulge)
	DNAbb = str(snakemake.params.dna_bulge)
	mm = str(snakemake.params.max_mismatch)
	PU = 'C'  # G = GPU C = CPU A = Accelerators -- I
	# === Wildcards ===
	guideref_name = str(snakemake.wildcards.guideref_name)
	fastaref_name = str(snakemake.wildcards.fastaref_name)

	# resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"

	# paths = listdir(resultsfolder)

	# Guide search params
	search_params = pickle.load(open(guide_search_params, 'rb'))
	# search_params = pickle.load(open(resultsfolder + "guide_search_params.pkl", 'rb'))

	# hg38 or consensus sequence
	fasta_fname = fasta_ref
	genome_name = fastaref_name
	# fasta_fname = '/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa'
	# genome_name = 'hg38'

	# hg38 guides found (but could be {alt_genome}_differences.csv
	guide_tab_fname = guides_report
	guides_src_name = guideref_name
	# guide_tab_fname = resultsfolder + 'Guides_found.csv'
	# guide_src_name = 'hg38'

	### Daniel---> Pycharm is not find subprocess.Popen(casoffinder...) without an absolute path. so I'm adding this
	# but I don't think its needed in the final version
	cas_off_expath = '/home/thudson/miniconda3/envs/edit/bin/cas-offinder'

	# defaults - we may allow users to change these cas-offinder settings?
	# according to Gorodkin et al. and Lin et al.  DNA bulges are even more tolerated than mismatches alone
	# https://www.nature.com/articles/s41467-022-30515-0
	# RNAbb = 0  # RNA bulge, a deletion in the off-target
	# DNAbb = 1  # DNA bulge, an insertion in the off-target
	# mm = 3  # max allowable mismatch
	# PU = 'C'  # G = GPU C = CPU A = Accelerators -- I don't really know which should be default?
	# casoff_params = (3, 0, 0, "C")
	casoff_params = (mm, RNAbb, DNAbb, PU)


if __name__ == "__main__":
	main()

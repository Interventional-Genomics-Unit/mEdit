import pickle
import pandas as pd
from subprocess import Popen
from os import listdir, remove
from collections import defaultdict

### This script takes ~45minutes to run and ~32seconds/guide found
# I recommend no more than 100 guides run at a time
'''
Cas-OFFinder 3.0.0 beta (Jan 24 2021)
Copyright (c) 2021 Jeongbin Park and Sangsu Bae
Website: http://github.com/snugel/cas-offinder

Usage: cas-offinder {input_filename|-} {C|G|A}[device_id(s)] {output_filename|-}
(C: using CPUs, G: using GPUs, A: using accelerators)
'''


def make_casoffinder_input(infile,fasta_fname,pam, pamISfirst, guidelen,guides,gnames,casoff_params):
    ## create input file for cas-offinder
    mm, RNAbb, DNAbb, PU = casoff_params

    with open(infile, 'w') as f:
        f.writelines(fasta_fname + "\n")
        line = 'N' * guidelen

        if pamISfirst:
            line = f"{pam}{line} {DNAbb} {RNAbb}" + "\n"
        else:
            line = f"{line}{pam} {DNAbb} {RNAbb}" + "\n"
        f.writelines(line)
        print(line)

        dpam = 'N' * len(pam)
        for grna, gname in zip(guides, gnames):
            if pamISfirst:
                line = f"{dpam}{grna} {mm} {gname}" + "\n"
            else:
                line = f"{grna}{dpam} {mm} {gname}" + "\n"
            f.writelines(line)
            print(line)


def cas_offinder_bulge(input_filename, output_filename,cas_off_expath,bulge):
    '''
     The cas-offinder off-line package contains a bug that doesn't allow bulges
    This script is partially a wrapper for cas-offinder to fix this bug
     created by...
    https://github.com/hyugel/cas-offinder-bulge

    '''
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
    outfn = fnhead + '_temp.txt'
    print("Running Cas-OFFinder (output file: %s)..." % outfn)
    p = Popen([cas_off_expath, casin, 'C', outfn])
    ret = p.wait()
    if ret != 0:
        print("Cas-OFFinder process was interrupted!")
        exit(ret)
    print("Processing output file...")

    with open(outfn) as fi, open(output_filename, 'w') as fo:
        fo.write('Guide_ID\tBulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n')
        for line in fi:
            entries = line.strip().split('\t')
            ncnt = 0

            if bulge == False:
                gid, mm = nobulge_dict[entries[0]]
                fo.write(f'{gid}\tX\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{entries[3]}\t{entries[4]}\t{entries[5]}\t0\n')
            else:
                if isreversed:
                    for c in entries[0][::-1]:
                        if c == 'N':
                            ncnt += 1
                        else:
                            break
                    if ncnt == 0:
                        ncnt = -len(entries[0])
                else:
                    for c in entries[0]:
                        if c == 'N':
                            ncnt += 1
                        else:
                            break

                if entries[0] in rnabulge_dic:
                    gid = id_dict[entries[0]]
                    for pos, query_mismatch, seq in rnabulge_dic[entries[0]]:
                        if isreversed:
                            tgt = (seq_pam + entries[0][len_pam:len_pam+pos] + seq + entries[0][len_pam+pos:-ncnt], entries[3][:len_pam+pos] + '-'*len(seq) + entries[3][len_pam+pos:-ncnt])
                        else:
                            tgt = (entries[0][ncnt:ncnt+pos] + seq + entries[0][ncnt+pos:-len_pam] + seq_pam, entries[3][ncnt:ncnt+pos] + '-'*len(seq) + entries[3][ncnt+pos:])
                        if query_mismatch >= int(entries[5]):
                            fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(gid,'RNA', tgt[0], tgt[1], entries[1], int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (isreversed and ncnt > 0 and entries[4] == "-") else 0), entries[4], int(entries[5]), len(seq)))
                else:
                    gid = id_dict[entries[0]]
                    nbulge = 0
                    if isreversed:
                        for c in entries[0][:-ncnt][len_pam:]:
                            if c == 'N':
                                nbulge += 1
                            elif nbulge != 0:
                                break
                        tgt = (seq_pam + entries[0][:-ncnt][len_pam:].replace('N', '-'), entries[3][:-ncnt])
                    else:
                        for c in entries[0][ncnt:][:-len_pam]:
                            if c == 'N':
                                nbulge += 1
                            elif nbulge != 0:
                                break
                        tgt = (entries[0][ncnt:][:-len_pam].replace('N', '-') + seq_pam, entries[3][ncnt:])
                    fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(gid,'X' if nbulge == 0 else 'DNA', tgt[0], tgt[1], entries[1], int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (isreversed and ncnt > 0 and entries[4] == "-") else 0), entries[4], int(entries[5]), nbulge))

        remove(fnhead + '_temp.txt')


def agg_results(output_filename,mmco):
    ots_dict = {}
    #grps =[[f'X_{n}'] for n in range(0,mmco)]
    with open(output_filename)as of:
        for line_out in of:
            entry = line_out.strip().split('\t')
            if entry[0] != 'Guide_ID':
                gid, btype, mm, bsize = entry[0], entry[1], entry[7], entry[8]
                if gid not in ots_dict.keys():
                    ots_dict[gid] ={}
                    for i in ['X','RNA','DNA']:
                        for j in range(mmco+1):
                            ots_dict[gid][(i,j)]= 0
                ots_dict[gid][(btype,int(mm))] += 1
    return ots_dict


def write_out_res(ots,gdf,casoff_params,resultsfolder,guide_tab_fname):
    df = pd.DataFrame(ots)

    #this is just a precautionary if this script was run more than once
    if 'Number of Mismatches' in gdf.columns:
        gdf = gdf.drop(columns='Number of Mismatches')

    #Update main guides table
    tot =df.sum(axis = 0)
    tot =tot.rename('Number of Mismatches')
    tot.index.name = 'Guide_ID'
    gdf = gdf.join(tot, on='Guide_ID')
    gdf.to_csv(guide_tab_fname, index = False)

    ## write out aggregate output
    df = df.reset_index()
    df = df.rename(columns = {'level_1': 'Number of Mismatches', 'level_0': 'BulgeType'})
    if casoff_params[1] == 0:
        df = df.loc[df.BulgeType != 'RNA']
    if casoff_params[2] == 0:
        df = df.loc[df.BulgeType != 'DNA']
    df = df.pivot_table(columns = ['BulgeType','Number of Mismatches'],aggfunc="sum")
    out = resultsfolder + 'Num_Mismatches.txt'
    df.to_csv(out, sep = '\t')


def run_casoffinder(resultsfolder,
                    fasta_fname,
                    guide_tab_fname,
                    search_params,
                    cas_off_expath,
                    genome_name,
                    guides_src_name,
                    casoff_params):
    gdf = pd.read_csv(guide_tab_fname)
    ots = {}
    gpr = gdf.groupby('Editor')
    if casoff_params[1:3] == (0, 0):
        bulge = False
    else:
        bulge = True
    # for each editor type find off_targets
    for editor, stats in gpr:
        infile = f"{resultsfolder}{genome_name}_{guides_src_name}_{editor}_casoffinder_input.txt"
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

        #start = time.time()
        ## with defualt setting 3mm, 1 dnabulge, 0 rnabulge,
        #run cas-offinder/ adjust input file for bulge
        output_filename = infile.replace('_input.txt', '_output.txt')
        cas_offinder_bulge(infile, output_filename, cas_off_expath, bulge)
        #end = time.time()
        #print(f'total_time for {editor}: {end - start}')
        #print(f'total_time per guide: {(end - start)/len(guides)}')

        #sum off-targets
        ot_dict = agg_results(output_filename,casoff_params[0])
        for k, v in ot_dict.items():
            ots[k] = v
    write_out_res(ots, gdf, casoff_params, resultsfolder, guide_tab_fname)


def main():
    # SNAKEMAKE IMPORTS
    # === Inputs ===
    guides_report = str(snakemake.input.guides_report_out)
    fasta_ref = str(snakemake.input.fasta_ref)
    guide_search_params = str(snakemake.input.guide_search_params)
    snv_site_info = str(snakemake.input.snv_site_info)
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

    resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"

    paths = listdir(resultsfolder)

    # Guide search params
    search_params = pickle.load(open(guide_search_params, 'rb'))

    # hg38 or consensus sequence
    fasta_fname = fasta_ref
    genome_name = fastaref_name

    # hg38 guides found (but could be {alt_genome}_differences.csv
    guide_tab_fname = guides_report
    guides_src_name = guideref_name

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
    casoff_params = (mm, RNAbb, DNAbb, PU)

    run_casoffinder(resultsfolder,
                    fasta_fname,
                    guide_tab_fname,
                    search_params,
                    cas_off_expath,
                    genome_name,
                    guides_src_name,
                    casoff_params)


if __name__ == "__main__":
    main()

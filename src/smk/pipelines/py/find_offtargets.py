import pickle
import pandas as pd
from subprocess import Popen, PIPE
from os import listdir, remove
from collections import defaultdict
from collections import OrderedDict
import scoring
from annotate import Transcript

###
###
#-----> DANIEL THIS SCRIPT ISNT DONE YET ********************* <-----------------------

### This script takes ~45minutes to run and ~32seconds/guide found
# I recommend no more than 100 guides run at a time
'''
Cas-OFFinder 3.0.0 beta (Jan 24 2021)
Copyright (c) 2021 Jeongbin Park and Sangsu Bae
Website: http://github.com/snugel/cas-offinder

Usage: cas-offinder {input_filename|-} {C|G|A}[device_id(s)] {output_filename|-}
(C: using CPUs, G: using GPUs, A: using accelerators)
'''

# INPUT LEG
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
        dpam = 'N' * len(pam)
        for grna, gname in zip(guides, gnames):
            if pamISfirst:
                line = f"{dpam}{grna} {mm} {gname}" + "\n"
            else:
                line = f"{grna}{dpam} {mm} {gname}" + "\n"
            f.writelines(line)



# INPUT/OUTPUT LEG
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


    p = Popen([cas_off_expath, casin, 'C', outfn])
    ret = p.wait()
    if ret != 0:
        print("Cas-OFFinder process was interrupted!")
        exit(ret)
    print("Processing output file...")

    # OUTPUT LEG
    with open(outfn) as fi, open(output_filename, 'w') as fo:
        fo.write('Coordinates\tDirection\tGuide_ID\tBulge type\tcrRNA\tDNA\tMismatches\tBulge Size\n')\
        #fo.write('Guide_ID\tBulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n')
        ot_coords = []
        for line in fi:
            entries = line.strip().split('\t')
            ncnt = 0


            if bulge == False:
                gid, mm = nobulge_dict[entries[0]]
                coord = f'{entries[1]}:{entries[2]}-{int(entries[2]) + len(entries[0])}'
                fo.write(f'{coord}\t{entries[4]}\t{gid}\tX\t{entries[0]}\t{entries[3]}\t{entries[5]}\t0\n')
                ot_coords.append(coord)
            else:
                if isreversed:
                    for c in entries[0][::-1]:
                        if c == 'N':
                            ncnt += 1
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

                            start = int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (isreversed and ncnt > 0 and entries[4] == "-") else 0)
                            coord = f'{entries[1]}:{start}-{int(start) + len(tgt[1])}'
                            ot_coords.append(coord)
                            fo.write(f'{coord}\t{entries[4]}\t{gid}\tRNA\t{tgt[0]}\t{tgt[1]}\t{int(entries[5])}\t{len(seq)}\n')

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
                    start =int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (isreversed and ncnt > 0 and entries[4] == "-") else 0)
                    btype = 'X' if nbulge == 0 else 'DNA'
                    coord = f'{entries[1]}:{start}-{int(start) + len(tgt[1])}'
                    ot_coords.append(coord)
                    fo.write(f'{entries[1]}:{start}-{start + len(tgt[1])}\t{entries[4]}\t{gid}\t{btype}\t{tgt[0]}\t{tgt[1]}\t{int(entries[5])}\t{nbulge}\n')
        remove(fnhead + '_temp.txt')
        remove(casin)
        editor = gid.split('_')[0]
        print(f'{len(ot_coords)} off targets found for {editor}')

def score_ot(crrna,otseq,editor):
    score = '.'
    if 'spCas9' in editor:
        # TODO: NEEDS
        raise Exception("The function scoring.cfd_score requires a path to the model files")
        score = scoring.cfd_score(crrna[:-3], otseq[:-4])
    return score



def annotate_ots(output_filename,annote_path):
    '''
    Reads output, Scores each off-target seq and annotates each off_target
    '''
    #annote_path = "/groups/clinical/projects/editability/tables/processed_tables/ncbiRefSeq.bed.gz"
    #output_filename = '/groups/clinical/projects/editability/medit_queries/medit_test/test_out/hg38_hg38_spCas9_casoffinder_output.txt'
    editor = output_filename.split('_casoff')[0].split('_')[-1]
    coords = []
    scores = []
    spec_scores = {}
    res = open(output_filename, 'r').readlines()
    for line in res[1:]:
        line = line.strip().split('\t')
        coords.append(line[0])
        if line[2] not in spec_scores.keys():
            spec_scores[line[2]] = 0
        score = score_ot(line[4], line[5],editor)
        if score != '.':
            if score != 1:
                spec_scores[line[2]] = spec_scores[line[2]] + score
        else:
            spec_scores[line[2]] = '.'
        scores.append(score)

    Transcript.load_transcripts(annote_path,coords)

    new_lines = []
    annotate_out = output_filename.replace('_output', '_annotated')
    with open(annotate_out, 'w') as anout:
        anout.write(res[0].strip() + f'\tAnnotation\tScore\n')
        cnt = 0
        for line in res[1:]:
            line = line.strip().split('\t')
            ans = Transcript.transcript(line[0])

            if ans != 'intergenic':
                tid, gname = ans.tx_info()[1:3]
                feature = ans.feature
                x = '|'.join([feature, gname, tid])
                new_line = line + [x, str(scores[cnt])]
                new_lines.append(new_line)

            else:
                x = 'intergenic'
                new_line = line + [x, str(scores[cnt])]
                new_lines.append(new_line)
            cnt += 1
            anout.write('\t'.join(new_line) + '\n')

    if editor == 'spCas9':
        for gid, sum_score in spec_scores.items():
            spec_scores[gid] = scoring.cfd_spec_score(sum_score)
    #remove(output_filename)
    return new_lines, spec_scores


def agg_results(lines,mmco):
    '''
    sums the number of off-targets for each guide
    aggregates by mismatch cutoff and bulge size
    '''
    ots_dict = {}
    for line in lines:
        gid, btype, mm, bsize = line[2], line[3], line[6], line[7]
        if gid not in ots_dict.keys():
            ots_dict[gid] ={}
            for i in ['X','RNA','DNA']:
                for j in range(mmco+1):
                    ots_dict[gid][(i,j)]= 0
        ots_dict[gid][(btype,int(mm))] += 1
    return ots_dict

# OUTPUT LEG
def write_out_res(ots,gdf,casoff_params,resultsfolder,guide_tab_fname):
    '''
    Adds two new columns to guides found table
    Creates and Aggregate summary text output
    '''
    df = pd.DataFrame(ots)#.rename(columns = {'0':'Off Target Score'})
    #this is just a precautionary if this script was run more than once
    if 'Num Off Targets per MM' in gdf.columns:
        gdf = gdf.drop(columns='Num Off Targets per MM')
    if 'Off Target Score' in gdf.columns:
        gdf = gdf.drop(columns='Off Target Score')

    # Update main guides table with specifity Scoress
    spec_scores = df.T.sort_index()['spec_score']
    df =df.drop(index = 'spec_score')

    #Add off_targets summary
    df.index = pd.MultiIndex.from_tuples(list(df.index), names=['BulgeType', 'Number of Mismatches'])
    df = df.reset_index()
    ot_per_mm = df.loc[df['BulgeType']=='X']
    ot_per_mm= ot_per_mm.iloc[:, 2:].T.sort_index()
    ot_per_mm=ot_per_mm.stack().astype('str').groupby(level=0).apply('|'.join)
    ot_per_mm=ot_per_mm.rename('Guide_ID')#('Num Off Targets per Mismatch')
    gdf = gdf.sort_values('Guide_ID')
    gdf['Off Target Score'] = list(spec_scores)
    gdf['Num Off Targets per MM'] = list(ot_per_mm)


    # create off_target summary of totals
    if casoff_params[1] == 0:
        df = df.loc[df.BulgeType != 'RNA']
    if casoff_params[2] == 0:
        df = df.loc[df.BulgeType != 'DNA']
    offtarget_summary_report = df.pivot_table(columns=['BulgeType', 'Number of Mismatches'], aggfunc="sum")

    #writeout
    sum_out = resultsfolder + 'OffTarget_Summary.txt'
    offtarget_summary_report.to_csv(sum_out, sep='\t')

    gdf.to_csv(guide_tab_fname, index = False)


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
        cas_offinder_bulge(infile, output_filename, cas_off_expath, bulge)

        new_lines, spec_scores = annotate_ots(output_filename,annote_path)

        ot_dict = agg_results(new_lines,casoff_params[0])#sum off-targets
        for k, v in ot_dict.items():
            v['spec_score'] = spec_scores[k]
            ots[k] = v
    write_out_res(ots, gdf, casoff_params, resultsfolder, guide_tab_fname)


def main():
    # SNAKEMAKE IMPORTS
    # === Inputs ===
    guides_report = str(snakemake.input.guides_report_out)
    fasta_ref = str(snakemake.input.fasta_ref)
    guide_search_params = str(snakemake.input.guide_search_params)
    snv_site_info = str(snakemake.input.snv_site_info)
    annote_path = str(snakemake.params.annote_path)
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

    #resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"

    paths = listdir(resultsfolder)

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
    #casoff_params = (3, 0, 0, "C")
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

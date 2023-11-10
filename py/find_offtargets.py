import pickle
import pandas as pd
from subprocess import Popen
from os import listdir, remove
from collections import defaultdict


'''
Cas-OFFinder 3.0.0 beta (Jan 24 2021)
Copyright (c) 2021 Jeongbin Park and Sangsu Bae
Website: http://github.com/snugel/cas-offinder

Usage: cas-offinder {input_filename|-} {C|G|A}[device_id(s)] {output_filename|-}
(C: using CPUs, G: using GPUs, A: using accelerators)
'''

#defaults
RNAbb = 1
DNAbb = 1
mm = 3
PU = 'C'
casoff_params = (mm,RNAbb,DNAbb,PU)

def make_casoffinder_input(infile,fasta_fname,pam, pamISfirst, guidelen,guides,gnames,casoff_params):
    mm, RNAbb, DNAbb, PU = casoff_params

    with open(infile, 'w') as f:
        f.writelines(fasta_fname + "\n")
        line = 'N' * guidelen

        if pamISfirst:
            line = f"{pam}{line} {RNAbb} {DNAbb}" + "\n"
        else:
            line = f"{line}{pam} {RNAbb} {DNAbb}" + "\n"
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
    print("Cas-OFFinder-bulge v1.2 beta (2016-07-06)")
    print("")
    print("Copyright (c) 2015 Jeongbin Park and Sangsu Bae")
    print("")
    print("Usage: cas-offinder-bulge {input_file} {C|G|A|D} {output_file}")
    print("(C: using CPUs, G: using GPUs, A: using accelerators, D: dry-run")
    fnhead = input_filename.replace("_input.txt", "")
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
            ids = []
            for raw_target, mismatch, gid in targets:
                print(raw_target, mismatch, gid)

                if isreversed:
                    target = raw_target.lstrip('N')
                    len_pam = len(raw_target) - len(target)
                    bg_tgts['N' * len_pam + target + 'N' * bulge_dna].add(mismatch)
                    for bulge_size in range(1, bulge_dna+1):
                        for i in range(1, len(target)):
                            bg_tgt = 'N' * len_pam + target[:i] + 'N' * bulge_size + target[i:] + 'N' * (bulge_dna - bulge_size)
                            bg_tgts[bg_tgt].add(mismatch)
                            ids.append(gid)

                    for bulge_size in range(1, bulge_rna+1):
                        for i in range(1, len(target)-bulge_size):
                            bg_tgt = 'N' * len_pam + target[:i] + target[i+bulge_size:] + 'N' * (bulge_dna + bulge_size)
                            bg_tgts[bg_tgt].add(mismatch)
                            rnabulge_dic[bg_tgt].append( (i, int(mismatch), target[i:i+bulge_size],gid) )
                            ids.append(gid)
                else:
                    target = raw_target.rstrip('N')
                    len_pam = len(raw_target) - len(target)
                    bg_tgts['N' * bulge_dna + target + 'N' * len_pam].add(mismatch)
                    for bulge_size in range(1, bulge_dna+1):
                        for i in range(1, len(target)):
                            bg_tgt = 'N' * (bulge_dna - bulge_size) + target[:i] + 'N' * bulge_size + target[i:] + 'N' * len_pam
                            bg_tgts[bg_tgt].add(mismatch)
                            ids.append(gid)

                    for bulge_size in range(1, bulge_rna+1):
                        for i in range(1, len(target)-bulge_size):
                            bg_tgt = 'N' * (bulge_dna + bulge_size) + target[:i] + target[i+bulge_size:] + 'N' * len_pam
                            bg_tgts[bg_tgt].add(mismatch)
                            rnabulge_dic[bg_tgt].append( (i, int(mismatch), target[i:i+bulge_size],gid) )
                            ids.append(gid)
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
                f.write(tgt + ' ' + str(max(mismatch)) + ' ' + ids[cnt] + '\n')
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
    p = Popen([f'{cas_off_expath}cas-offinder', casin, 'C', outfn])
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
                    for pos, query_mismatch, seq, gid in rnabulge_dic[entries[0]]:
                        if isreversed:
                            tgt = (seq_pam + entries[0][len_pam:len_pam+pos] + seq + entries[0][len_pam+pos:-ncnt], entries[3][:len_pam+pos] + '-'*len(seq) + entries[3][len_pam+pos:-ncnt])
                        else:
                            tgt = (entries[0][ncnt:ncnt+pos] + seq + entries[0][ncnt+pos:-len_pam] + seq_pam, entries[3][ncnt:ncnt+pos] + '-'*len(seq) + entries[3][ncnt+pos:])
                        if query_mismatch >= int(entries[5]):
                            fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(gid,'RNA', tgt[0], tgt[1], entries[1], int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (isreversed and ncnt > 0 and entries[4] == "-") else 0), entries[4], int(entries[5]), len(seq)))
                else:
                    bulge = 0
                    if isreversed:
                        for c in entries[0][:-ncnt][len_pam:]:
                            if c == 'N':
                                bulge += 1
                            elif bulge != 0:
                                break
                        tgt = (seq_pam + entries[0][:-ncnt][len_pam:].replace('N', '-'), entries[3][:-ncnt])
                    else:
                        for c in entries[0][ncnt:][:-len_pam]:
                            if c == 'N':
                                bulge += 1
                            elif bulge != 0:
                                break
                        tgt = (entries[0][ncnt:][:-len_pam].replace('N', '-') + seq_pam, entries[3][ncnt:])
                    fo.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(gid,'X' if bulge == 0 else 'DNA', tgt[0], tgt[1], entries[1], int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (isreversed and ncnt > 0 and entries[4] == "-") else 0), entries[4], int(entries[5]), bulge))

        remove(fnhead + '_temp.txt')



def agg_results(output_filename):

    ots_dict = {}
    '''
    with open(input_filename)as inf:
        ots_sub = {}
        for line_in in inf:
            entry = line_in.strip().split(' ')
            if len(entry) > 2 and len(entry[-1]) < 3:
                RNAbb,DNAbb = int(entry[-2]),  int(entry[-1])
            if len(entry) > 2 and len(entry[-1]) > 3:
                x = [0] * int(entry[1])
                ots_sub['X_0'] = x
                for i in range(1,RNAbb+1):
                    ots_sub[f'RNA_{i}'] = x
                for i in range(1,DNAbb+1):
                    ots_sub[f'DNA_{i}'] = x
                ots_dict[entry[-1]] = ots_sub
        '''

    with open(output_filename)as of:
        for line_out in of:
            entry = line_out.strip().split('\t')
            if entry[0] != 'Guide_ID':
                gid, btype, mm, bsize = entry[0], entry[1], entry[7], entry[8]
                if gid not in ots_dict.keys():
                    ots_dict[gid] = {'X':0,'RNA':0,'DNA':0}
                ots_dict[gid][btype] += 1
    return ots_dict


def run_casoffinder(gdf, fasta_fname,resultsfolder,genome_name,cas_off_expath,casoff_params):
    infiles = []
    ots = {}
    gpr = gdf.groupby('Editor')
    if casoff_params[1:3] == (0, 0):
        bulge = False
    else:
        bulge = True

    for editor,stats in gpr:
        infile = f"{resultsfolder}{genome_name}{editor}_casoffinder_input.txt"
        infiles.append(infile)
        pam, pamISfirst, guidelen = search_params[editor][0:3]
        guides, gnames = list(stats.gRNA), list(stats.Guide_ID)

        # make input file
        make_casoffinder_input(infile, fasta_fname,pam, pamISfirst, guidelen, guides, gnames,casoff_params)

        #run cas-offinder/ adjust input file for bulge
        output_filename = infile.replace('_input.txt', '_output.txt')
        cas_offinder_bulge(infile, output_filename, cas_off_expath,bulge)

        #sum off-targets
        ot_dict = agg_results(output_filename)
        for k, v in ot_dict.items():
            ots[k] = v
    gdf.join(pd.DataFrame(ots).T, on = 'Guide_ID')

def get_offtargets(resultsfolder,fasta_fname,hg38guide_results,refgenome_name,altguide_results,searchp_path):
    hg38_gdf = pd.read_csv(hg38guide_results)
    if altguide_results != False:





resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"

##hg38 or consensus sequence
fasta_fname = "/groups/clinical/projects/medit_analysis/private/consensus_refs/hg38/PG_WGS_HG38.fa"
genome_name = altguide_results[0].split('/')[-1].split('.')[0]

# hg38 guides found
hg38guide_results = '/groups/clinical/projects/editability/medit_queries/medit_test/test_out/2023-11-01_Guides_found.csv'
refgenome_name = 'HG38'

# Alt guides found from process_genome.py
altguide_results =[x for x in paths if x.endswith('differences.csv')]
if len(altguide_results) == 0:
    altguide_results = False
    altgenome_name = False
else:
    altguide_results  = altguide_results[0]
    altgenome_name = altguide_results[0].split('/')[-1].split('_')[0]

searchp_path = [resultsfolder + x for x in paths if x.endswith('guide_search_params')][0] #Get guide search params used
search_params = pickle.load(open(searchp_path, 'rb'))



### Daniel---> Pycharm is not find subprocess.Popen(casoffinder...) without an absolute path. so I'm adding this
# but I don't think its needed in the final version
cas_off_expath = '/home/thudson/miniconda3/envs/edit/bin/'


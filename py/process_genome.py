import pandas as pd
import gzip
from Bio import SeqIO
#################################
# takes a pangenome GTF file and maps to hg38 refseq annotations by EnsembleIDs
#################################

out_folder = '/groups/clinical/projects/editability/tables/processed_tables/'
hg38_refseq_path = "/groups/clinical/projects/editability/tables/processed_tables/ncbiRefSeq.txt.gz"
gtf_path = "/groups/clinical/projects/editability/genomes/pangenomes/HG02257/paternal/GCA_018466835.1.gtf.gz"


def get_eids(ch,hg38_refseq_path):
    '''
    extract ensemble ids in refernece genome(hg38)
    '''
    eids = []
    lines = []
    for line in gzip.open(hg38_refseq_path, 'rt'):
        tokens = line.split('\t')
        if tokens[2] == f'chr{ch}':
            if tokens[0] != '-':
                lines.append(tokens)
                eids.append(tokens[0])
    return eids, lines


def extract_eids(eids, gtf_path):
    '''
     extract match ensemble ids in pangenome gtf file to those found in refernece genome(hg38)
    '''
    names = ['Chr', 'source', 'feature', 'Start', 'End', 'score', 'strand', 'frame', 'attribute']
    lines = []
    chrs = []
    for line in gzip.open(gtf_path, 'rt'):
        if line.startswith('#') == False:
            tokens = line.split("\t")
            entry = dict(zip(names, tokens))
            chrs.append(entry['Chr'])
            if entry['feature'] != 'gene':
                a = entry['attribute']
                e = a[a.find("ENST0000"):a.find("ENST0000") + 17]
                if e in eids:
                    lines.append(tokens)
    return lines


def create_txt_file(ref_lines, gtf_lines, ch):
    '''
     create a formated refseq txt file for pangenome
    '''
    labels = ['eid', 'tid', 'chrom', 'strand', 'txStart', 'txEnd',
              'cdsStart', 'cdsEnd', 'exonStarts', 'exonEnds',
              'name', 'exonFrames', '5utrstarts', '5utrends', '3utrstarts', '3utrends']

    attr_names = ['projection_parent_transcript', 'exon_number', 'gene_name']
    new_txts = {}
    for x in gtf_lines:
        attribute = x[-1].split(';')
        found = ['-'] * len(attr_names)
        for n in attr_names:
            for a in attribute:
                if a.strip().split(" ")[0] == n:
                    found[attr_names.index(n)] = a.strip().split(" ")[1].replace("\"", "")

        if x[2] == 'transcript':
            new_txts[found[0]] = [found[0], '-', f'chr{ch}|{x[0]}', x[6], x[3], x[4], '-', '-', '-', '-', found[2], '-',
                                  '-', '-', '-', '-']

        if found[0] in new_txts.keys():
            if x[2] == 'exon':
                new_txts[found[0]][8] = new_txts[found[0]][8] + ',' + x[3]
                new_txts[found[0]][9] = new_txts[found[0]][9] + ',' + x[4]
            elif x[2] == 'CDS':
                new_txts[found[0]][6] = x[3]
                new_txts[found[0]][7] = x[4]
            elif x[2] == 'five_prime_utr':
                new_txts[found[0]][12] = x[3]
                new_txts[found[0]][13] = x[4]
            elif x[2] == 'three_prime_utr':
                new_txts[found[0]][14] = x[3]
                new_txts[found[0]][15] = x[4]
            else:
                pass
    df = pd.DataFrame(new_txts.values(), columns=labels)
    df.exonStarts = df.exonStarts.str.replace('-,', "")
    df.exonEnds = df.exonEnds.str.replace('-,', "")
    df = df.drop(columns=['exonFrames', 'tid'])
    refdf = pd.DataFrame(ref_lines, columns=labels[0:-4])
    refdf = refdf[['eid', 'tid', 'exonFrames']]
    df = df.join(refdf.set_index('eid'), on='eid')
    df = df[labels[0:-4]]
    return df


def gtf_to_refseq(gtf_path,hg38_refseq_path,out_folder):

    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
              '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'Y', 'X']
    df = pd.DataFrame()

    for ch in chroms:
        print(ch)
        eids, ref_lines = get_eids(ch,hg38_refseq_path)
        gtf_lines = extract_eids(eids, gtf_path)
        temp = create_txt_file(ref_lines, gtf_lines, ch)
        df = pd.concat([df,temp])
    x = df.to_dict('tight')['data']
    out = out_folder + gtf_path.split('/','')[0].replace('.gtf.gz','txt.gz')
    outfile = gzip.open(out,'wt')
    for line in x:
        outfile.write('\t'.join(line))
    outfile.close()

###test
def get_refseq_entry(term, field,annote_path,fasta_path):
    '''
    Using ncbiRefSeq.txt to find cds features by either interval, gene name or transcript ID
    example input:
    term, field = 'NM_000532.5', 'tid'
    term, field = 'ENST00000251654.9', 'eid'
    term, field = 'PCCB','name'
    term,field =  'chr3:136327650-136327650','interval'
    '''

    labels = ['eid', 'tid', 'chrom', 'strand', 'txStart', 'txEnd',
              'cdsStart', 'cdsEnd', 'exonStarts', 'exonEnds', 'name', 'exonFrames']


    if field != 'interval':
        not_found = True
        for line in gzip.open(annote_path, 'rt'):
            tokens = line.split('\t')
            entry = dict(zip(labels, tokens))
            if term in entry[field]:
                not_found = False
                break

    else:  # only used for intervals search
        not_found = True
        ch = term.split(":")[0]
        start, end = term.split(":")[1].split('-')
        pos = int((int(start) + int(end)) / 2)

        for line in gzip.open(annote_path, 'rt'):
            tokens = line.split('\t')
            entry = dict(zip(labels, tokens))
            if ch in entry['chrom']:
                if pos in range(int(entry['txStart']), int(entry['txEnd'])):
                    not_found = False
                    break

    if not_found:
        entry = None
        print(f"{term} not found in refseq data")
    else:
        fasta_seq = SeqIO.parse(gzip.open(fasta_path, 'rt'), 'fasta')
        for record in fasta_seq:
            if entry['chrom'].split('|')[1] in record.id:
                print(record.id)
                fasta_seq = record.seq
                break
        tx_seq = fasta_seq[int(entry['txStart']):int(entry['txEnd'])]
    return entry, tx_seq

term,field =  'NM_000532.5', 'tid'
annote_path = '/groups/clinical/projects/editability/tables/processed_tables/HG02257.paternal.f1.txt.gz'
fasta_path = '/groups/clinical/projects/editability/genomes/test_pangenome/HG02257.paternal.f1_assembly_v2_genbank.fa.gz'

entry,tx_seq = get_refseq_entry(term,field,annote_path,fasta_path)


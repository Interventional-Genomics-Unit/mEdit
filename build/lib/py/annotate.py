import gzip
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq


def get_refseq_entry(term, field, annote_path):
    '''
    Using ncbiRefSeq.txt to find cds features by either interval, gene name or transcript ID
    example input:
    term, field = 'NM_000532.5', 'tid'
    term, field = 'ENST00000251654.9', 'eid'
    term, field = 'PCCB','name'
    term,field =  'chr3:136250339-136330169','interval'
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

        if not_found:
            entry = None
            print(f"{term} not found in refseq data")

            if '.' in term:
                new_term = term.split('.')[0]
                print(f'searching for {new_term} instead')
                entry = get_refseq_entry(new_term, field,annote_path)

    else:  # only used for intervals search
        not_found = True
        ch = term.split(":")[0]
        start, end = term.split(":")[1].split('-')
        pos = int((int(start) + int(end)) / 2)

        for line in gzip.open(annote_path, 'rt'):
            tokens = line.split('\t')
            entry = dict(zip(labels, tokens))
            if ch == entry['chrom']:
                if pos in range(int(entry['txStart']), int(entry['txEnd'])):
                    not_found = False
                    break
        if not_found:
            entry = None
            print(f"{term} not found in refseq data")

    return entry


def get_cds_info(tx_seq, entry):
    '''
    uses entry info to find cds (without utr's)
    '''
    exon_starts = entry['exonStarts'][:-1].split(',')
    exon_ends = entry['exonEnds'][:-1].split(',')
    exon_frames = entry['exonFrames'].replace("\n", "")[:-1].split(',')
    tx_start = int(entry['txStart'])

    exons = [(int(exon_starts[i]) - tx_start, int(exon_ends[i]) - tx_start) for i in range(len(exon_ends))]
    for i in range(len(exon_frames)):
        if exon_frames[i] == '-1':  # -1 means entire exon is UTR
            exons = exons[1:]
            exon_starts = exon_starts[1:]
        else:
            break
    for i in range(1, len(exon_frames)):
        if exon_frames[-i] == '-1':
            exons = exons[0:len(exons) - 1]
            exon_ends = exon_ends[0:-1]
        else:
            break

    # Determine the stop and start of UTR
    if len(exons) > 0:
        exons[0] = (int(entry['cdsStart']) - int(exon_starts[0]) + exons[0][0], exons[0][1])
        exons[-1] = (exons[-1][0], exons[-1][1] - (int(exon_ends[-1]) - int(entry['cdsEnd'])))

        cds = Seq(''.join([str(tx_seq)[a:b] for a, b in exons]))
        if entry['strand'] == '-':
            cds = cds.reverse_complement()
    else:
        cds = None
        exons = None

    # translation = cds.translate()
    return [exons, tx_seq, cds]

def find_transcript_info(term, fasta, annote_path):
    '''
    Using a Refseq Transcript_ID, Ensembl Transcript_ID or coordinates find transcript annotations and transcript sequence
    from either a genome fasta path or given genome sequence
    '''
    # id= 'NM_000532.5' or 'ENST00000251654.9'
    # fasta = f"/groups/clinical/projects/clinical_shared_data/hg38/hg38_chr20.fa.gz"
    if type(fasta) == str:
        fasta_seq = SeqIO.read(gzip.open(fasta, 'rt'), 'fasta')
    else:
        fasta_seq = fasta
    field = 'eid' if term.startswith('E') else 'tid' if term.startswith('N') else 'interval'

    entry = get_refseq_entry(term=term, field=field, annote_path=annote_path)
    if entry != None:
        tx_seq = fasta_seq.seq[int(entry['txStart']):int(entry['txEnd'])]
        tid_info = get_cds_info(tx_seq, entry)
    else:
        tid_info = None
    return entry, tid_info


def find_codons(dist_from_cds_start):
    '''
    Finds reading frame of SNV in extracted sequence
    '''
    rf = 1 if dist_from_cds_start % 3 == 2 else 2 if dist_from_cds_start % 3 == 0 else 0
    return rf


def find_snvseq_info(seq,snvpos,exons, cds, entry,window):
    # returns - sequence,feature,translation(if needed)
    # feature: non-coding, utr5,ut3,intron,exon, start_codon, stop_codon
    # snvpos, alt = 11576257, 'T'
    global dist_from_cds_start
    feature, rf = None, None

    t_snvpos = int(snvpos) - int(entry['txStart'])
    cdstart, cdsend = int(entry['cdsStart']) - int(entry['txStart']), int(entry['txEnd']) - int(entry['txStart'])
    strand = entry['strand']

    if t_snvpos < 0 or cds == None:
        # not in transcript - shouldn't happen or else no entry would be found
        feature = 'non-coding'

    else:
        if t_snvpos in range(cdstart, cdsend + 1):
            # in CDS
            # find if utr
            feature = 'intron'
            if t_snvpos in range(cdstart, exons[0][0] + 1):
                feature = '3utr' if strand == '-' else '5utr'
            elif t_snvpos in range(exons[-1][1], cdsend + 1):
                feature = '5utr' if strand == '-' else '3utr'

            # find if exon or intron
            else:
                exon_n = 0
                for x in exons:
                    # if in exon find reading frame
                    if t_snvpos in range(x[0], x[1] + 1):
                        # stop and start codon
                        feature = 'exon'
                        dist = sum([e[1] - e[0] for e in exons[0:exon_n]])
                        dist_from_cds_start = dist + (t_snvpos - x[0])

                        if strand == '-':
                            dist_from_cds_start = (len(cds) - dist_from_cds_start) + 1

                        if dist_from_cds_start < 3:
                            feature = 'start_codon'

                        if dist_from_cds_start > len(cds) - 3:
                            feature = 'stop_codon'

                        rf = find_codons(dist_from_cds_start)
                        break
                    exon_n += 1
        else:
            # in transcript but not in cds
            seq = str(seq)

            if len(SeqUtils.nt_search(seq[window - 6:window + 5], 'TTTATT')) > 1 or len(
                    SeqUtils.nt_search(seq[window - 6:window + 5], 'AATAAA')) > 1:
                feature = 'polya'
            elif len(SeqUtils.nt_search(seq[window - 6:window + 5], 'TATAAA')) > 1 or len(
                    SeqUtils.nt_search(seq[window - 6:window + 5], 'ATATTT')) > 1:
                feature = 'promoter'
            elif len(SeqUtils.nt_search(seq[window - 7:], 'GGNCAATCT')) > 1:
                if len(SeqUtils.nt_search(seq[window - 7:window + 6], 'GGNCAATCT')) > 1:
                    feature = 'promoter'
                else:
                    feature = 'TSS'
            elif len(SeqUtils.nt_search(seq[window + 6:], 'AGATTGNCC')) > 1:
                if len(SeqUtils.nt_search(seq[window - 7: window + 6], 'AGATTGNCC')) > 1:
                    feature = 'promoter'
                else:
                    feature = 'TSS'
            else:
                feature = 'flanking'

    return feature, rf


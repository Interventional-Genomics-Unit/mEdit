import pandas as pd
import re
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SeqUtils
import math


# tables and dictionaries
def check_queries(queries):
    '''
    Extracts Prefix and suffix and Validates HGVS terms given by user
    returns dict of prefix and suffix of valdiated terms
    ex. prefix =

    '''
    # verify query syntax
    rsuffix = r"[0-9cm.]*[\d\\-\\+\\*]*[\w>\\+\\-\\*]*>\w+"
    rprefix = r"((N(M|G|C|R)_\d+)|(m))."
    new_queries, prefixes, suffixes = [], [], []
    for q in set(queries):
        if re.search(rsuffix, q) :
            if re.search(rprefix,q):
                prefixes.append(re.search(rprefix, q).groups()[0])
                suffixes.append(re.search(rsuffix, q).group(0))
                # suffixes.append(re.search(rsuffix, q).captures()[0])
                new_queries.append(q)
    n = len(new_queries)
    print(f'{n} out of {len(queries)} validated')
    if n == 0:
        print('Query are not in the correct HGVS Format')
    return new_queries, prefixes, suffixes


def get_chroms_from_RefID(processed_tables,queries):
    # find chrom from HGVS Lookup table
    HGVSlookup_path = f"{processed_tables}/HGVSlookup.csv"
    newchromsd = {}

    hgvs_tab = pd.read_csv(HGVSlookup_path)
    queries, prefixes, suffixes = check_queries(queries)
    matches = hgvs_tab.loc[hgvs_tab['HGVS'].isin(prefixes)]
    chromsd = pd.DataFrame(data = list(zip(queries, prefixes, suffixes)),
                            columns = ['queries','HGVS','suffixes']).join(matches.set_index('HGVS'),on = 'HGVS')
    chromsd = chromsd.drop_duplicates(subset= ['queries'])
    chromsd = chromsd.to_dict('split')['data']

    # add hgvs prefix and suffix dictionary with chroms in keys
    for d in chromsd:  # chomsd =  'data': [['NM_000518', '11'],...
        if d[3] not in newchromsd.keys():
            newchromsd[d[3]] = [(d[1], d[2])]
        else:
            newchromsd[d[3]].append((d[1],d[2]))
    return newchromsd


def get_seqinfo(queries,qtype, datadir):

    processed_tables = f"{datadir}/processed_tables"
    guidelookup_tables = f"{datadir}/processed_tables/guide_acquisition_tables"
    chromsd = get_chroms_from_RefID(processed_tables,queries)
    hgvs_info = pd.DataFrame()

    for chrom, terms in chromsd.items():
        tempdf = pd.read_csv(f'{guidelookup_tables}/{chrom}_guide_variant.txt')
        hgvs = []
        for t in terms:
            prefix, suffix= t[0], t[1]
            ids = tempdf.loc[tempdf['HGVS_ID'].str.startswith(prefix), 'HGVS_ID']
            hgvs.append([h for h in ids if suffix in h][0])
        tempdf = tempdf.loc[tempdf['HGVS_ID'].isin(hgvs)]
        if len(tempdf['HGVS_ID']) > 0:
            hgvs_info = pd.concat([tempdf, hgvs_info])
    variantseq_dict = hgvs_info.set_index('HGVS_ID').to_dict('index')


    return variantseq_dict
    #NM_000518.5(HBB):c.114G>A (p.Trp38Ter) {'AlleleID': 30444, 'RefAlleleVCF': 'C', 'AltAlleleVCF': 'T', 'GeneID': 3043,
    # 'GeneSymbol': 'HBB', 'Chr': 11, 'PositionVCF': 5226778, 'Start': 5225464.0, 'End': 5227071.0, 'Strand': '-',
    # 'ChrID': 'NC_000011.10', 'TranscriptID': 'NM_000518.5', 'ProteinID': 'NP_000509.1', 'cdsStart': 5225597.0,
    # 'cdsEnd': 5227021.0, 'exonStart': '5225463,5226576,5226929,', 'MC': 'nonsense', 'Feature': 'exon'}


class DataHandler():
    '''
    search fasta for guides given reqs
    pam, len, targe table window and optionally scores
    '''

    def __init__(self, hgvs_id: str, data: dict, fasta_path):
        '''
        :param hgvs_id: ex: 'NM_000532.5(PCCB):c.1316A>G (p.Tyr439Cys)'
        :param data: ex: {'AlleleID': 227253, 'RefAlleleVCF': 'A', 'AltAlleleVCF': 'G', 'GeneID': 5096, 'GeneSymbol': 'PCCB', 'Chr': 3, 'Strand': '+',
        'PositionVCF': 136327650, 'TranscriptID': 'NM_000532.5', 'cdsStart': 136250375.0, 'cdsEnd': 136330026.0,
        'exonStart': '136250339,136255855,136256554,136260478,136261951,136283836,136293755,136297951,136301029,136316940,136326802,136327154,136327633,136328757,136329904,',
         'MC': 'missense_variant', 'Feature': 'exon'}

        :param fasta_path:
        '''
        #paths
        self.fasta_path = fasta_path

        #search data
        self.hgvs_id= hgvs_id
        self.chrom = data['Chr']
        self.NC_ref_allele = data['RefAlleleVCF']
        self.NC_alt_allele = data['AltAlleleVCF']
        self.strand = data['Strand']
        self.NM_ref_allele = data['RefAlleleVCF'] if data['Strand'] == '+' else str(Seq(data['RefAlleleVCF']).complement())
        self.NM_alt_allele = data['AltAlleleVCF'] if data['Strand'] == '+' else str(Seq(data['AltAlleleVCF']).complement())
        self.SNV_chr_pos = data['PositionVCF']
        self.alleleid = data['AlleleID']
        self.mc = data['MC']

        #search params
        self.pam = str() # Ex. 'NGG'
        self.pamISfirst = False #Boolean
        self.win_size = list() # Ex. list [4,8]
        self.scoring = None #Ex. True/False
        self.guidelen = 20

        #variables
        self.search_window = 30
        self.extracted_seq = str()

        #outputs
        self.guides_found ={'HGVS_ID':[],'Editor': [],'Guide_ID':[],'Chr': [],
                            'Start': [],'End': [],
                            'Strand': [],'gRNA': [],
                            'Pam': [],'Score':[],
                            'SNV Position': [], 'Variant_Molecular_Consequence':[]}

        self.BEguides_found ={'HGVS_ID': [],
                     'Base Editor': [],'Guide_ID': [],'Chr': [],
                     'Start': [],
                     'End': [],
                     'gRNA': [],
                     'Pam': [],
                     'SNV Position': [],
                     'Codon Position': [],
                     'Strand': [],
                     'Reference (Codon>AA)': [],
                     'Alternate (Codon>AA)': [],
                     'BE Converted (Codon>AA)': [],
                     'Conversion Type': [],
                     'Bystander': [],
                     'Variant_Molecular_Consequence': []
                     }

        #tables and dictionaries
        self.doenchParams = [
            # pasted/typed table from PDF and converted to zero-based positions
            (1, 'G', -0.2753771), (2, 'A', -0.3238875), (2, 'C', 0.17212887), (3, 'C', -0.1006662),
            (4, 'C', -0.2018029), (4, 'G', 0.24595663), (5, 'A', 0.03644004), (5, 'C', 0.09837684),
            (6, 'C', -0.7411813), (6, 'G', -0.3932644), (11, 'A', -0.466099), (14, 'A', 0.08537695),
            (14, 'C', -0.013814), (15, 'A', 0.27262051), (15, 'C', -0.1190226), (15, 'T', -0.2859442),
            (16, 'A', 0.09745459), (16, 'G', -0.1755462), (17, 'C', -0.3457955), (17, 'G', -0.6780964),
            (18, 'A', 0.22508903), (18, 'C', -0.5077941), (19, 'G', -0.4173736), (19, 'T', -0.054307),
            (20, 'G', 0.37989937), (20, 'T', -0.0907126), (21, 'C', 0.05782332), (21, 'T', -0.5305673),
            (22, 'T', -0.8770074), (23, 'C', -0.8762358), (23, 'G', 0.27891626), (23, 'T', -0.4031022),
            (24, 'A', -0.0773007), (24, 'C', 0.28793562), (24, 'T', -0.2216372), (27, 'G', -0.6890167),
            (27, 'T', 0.11787758), (28, 'C', -0.1604453), (29, 'G', 0.38634258), (1, 'GT', -0.6257787),
            (4, 'GC', 0.30004332), (5, 'AA', -0.8348362), (5, 'TA', 0.76062777), (6, 'GG', -0.4908167),
            (11, 'GG', -1.5169074), (11, 'TA', 0.7092612), (11, 'TC', 0.49629861), (11, 'TT', -0.5868739),
            (12, 'GG', -0.3345637), (13, 'GA', 0.76384993), (13, 'GC', -0.5370252), (16, 'TG', -0.7981461),
            (18, 'GG', -0.6668087), (18, 'TC', 0.35318325), (19, 'CC', 0.74807209), (19, 'TG', -0.3672668),
            (20, 'AC', 0.56820913), (20, 'CG', 0.32907207), (20, 'GA', -0.8364568), (20, 'GG', -0.7822076),
            (21, 'TC', -1.029693), (22, 'CG', 0.85619782), (22, 'CT', -0.4632077), (23, 'AA', -0.5794924),
            (23, 'AG', 0.64907554), (24, 'AG', -0.0773007), (24, 'CG', 0.28793562), (24, 'TG', -0.2216372),
            (26, 'GT', 0.11787758), (28, 'GG', -0.69774)]

    def set_guide_search_params(self,pam, pamISfirst,win_size,scoring,guidelen):
        self.pam = pam
        self.pamISfirst = pamISfirst
        self.win_size = win_size
        self.scoring = scoring
        self.guidelen = guidelen

    def calcDoenchScores(self, seq):
        """
        Input is a 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'
        """
        intercept = 0.59763615
        gcHigh = -0.1665878
        gcLow = -0.2026259

        assert (len(seq) == 30)
        score = intercept
        guideSeq = seq[4:24]
        gcCount = guideSeq.count("G") + guideSeq.count("C")

        if gcCount <= 10:
            gcWeight = gcLow
        if gcCount > 10:
            gcWeight = gcHigh

        score += abs(10 - gcCount) * gcWeight

        for pos, modelSeq, weight in self.doenchParams:
            subSeq = seq[pos:pos + len(modelSeq)]
            if subSeq == modelSeq:
                score += weight
        score = int(100 * (1.0 / (1.0 + math.exp(-score))))

        return score
    '''
    def cas_offinder(input, output):
        cmdArgs2 = ['cas-offinder', input, 'C', output]
        print(cmdArgs2)
        call(cmdArgs2)
        # output.flush()
        columns = ['crRNA', 'Chromosome', 'Position', 'Off-targets', 'Strand', 'No of Mismatches']
        ots = pd.read_csv(output, sep="\t", header=None, names=columns)
        # ots
        ots['Chromosome'] = ots['Chromosome'].str.split().str[0]

        ots = ots.sort_values(by='No of Mismatches')
        return (ots)
    '''

    def get_AAconversion_type(self,codon1,codon2):
        '''
        codon1: codon of Alt allele to be changed by BE
        codon2: codon after changed by BE
        '''
        aa_groups = [["G","A","V","L","I","M","F","Y","W"],
                     ["S","Q","T","N"],
                     ["C","G","P"],
                     ["D","E"],
                     ["K","H","R","Q"]]
        codon1,codon2 = Seq(codon1),Seq(codon2)
        aa1, aa2 = codon1.translate(), codon2.translate()
        mtype = ""
        if aa1 == aa2:
            mtype = 'Synonymous'
            if codon1 == codon2:
                mtype= 'Silent'
        else:
            if codon2 in ['TAA', 'TAG', 'TGA']:
                mtype = 'Nonsense'
            elif [aa_groups.index(x) for x in aa_groups if str(aa1) in x] == [aa_groups.index(x) for x in aa_groups if str(aa2) in x]:
                mtype = 'Conservative'
            else:
                mtype = 'Non-conservative'
        return mtype

    def codingType(self):
        '''
        Determine if translated region
        '''

        #genenomic_strand = self.vardf.strand.iloc[0]
        exon_pos = re.search("([\d\\-\\+]*)\w+>\w+", self.hgvs_id).groups()[0] #Ex. 1541+1 or 141

        if "-" in exon_pos or "+" in exon_pos:
            if exon_pos.startswith("-") or exon_pos.startswith("+"):
                coding_type = "5UTR" if exon_pos.startswith("-") else "3UTR"
            else:
                coding_type = 'intron'
        else:
            coding_type = "exon"

        return coding_type, exon_pos

    def getBE(self,guides,conversion,win_size,name):
        '''
        Finds codon level SNV and determines if the Base Editor Conversion can work
        *TO DO : Check HGVS term to find protein change and then bypass this 'manual' conversion

        '''
        f_delta = conversion # CT or AG
        r_delta = [str(Seq(conversion[0]).complement()),str(Seq(conversion[1]).complement())] ## the complment so GA or TC

        #Determine Reading frame
        coding_type, exon_pos = self.codingType()
        coding_strand = self.strand
        ref_allele = self.NM_ref_allele
        alt_allele = self.NM_alt_allele

        n =1
        for i in range(len(guides)):
            editor, guide, pam, strand, snvpos, score, start, end = guides[i]

            target_bases = Seq(guide[win_size[0] - 1:win_size[1] + 1])  # Bases inside 4-8 window

            # Converted case
            convert = str(f_delta[1] if alt_allele == f_delta[0] else r_delta[1])
            bystander = target_bases.count(f_delta[0]) - 1

            if coding_type != 'exon':
                alt = alt_allele
                ref = ref_allele
                codonpos = 'NA'
                ctype = 'NA'
                self.add_BEguides(name, guide, pam, strand, snvpos, score, start, end,codonpos,ref,alt,convert,ctype,bystander,n)
                n += 1

            else:
                ## In Exon
                #Determine codon and tranlated product
                coding_pos = round(int(exon_pos) /3,1) #find exon positiom
                if strand == coding_strand:
                    rf = 2 if coding_pos - round(coding_pos) == 0 else 1 if (coding_pos - round(
                        coding_pos)) == 0.5 else 0
                else:
                    rf = 0 if coding_pos - round(coding_pos) == 0 else 1 if (coding_pos - round(
                        coding_pos)) == 0.5 else 2

                codon_start = int(int(snvpos) - rf)

                #Alternative codon
                alt_codon = Seq(guide[codon_start:codon_start+3])
                if coding_strand != strand:
                    alt_codon = alt_codon.reverse_complement()
                aa_alt = alt_codon.translate()

                #Reference codon
                ref_codon = Seq("".join(alt_codon[x] if x != rf else ref_allele for x in [0,1,2]))
                if coding_strand != strand:
                    ref_codon = Seq("".join(alt_codon[x] if 2-x != rf else ref_allele for x in [0, 1, 2]))
                aa_ref = ref_codon.translate()

                ##Converted allele
                convert = f_delta[1] if alt_allele == f_delta[0] else r_delta[1]
                new_codon = Seq("".join(alt_codon[x] if x != rf else convert for x in [0,1,2]))
                if coding_strand != strand:
                    new_codon = Seq("".join(alt_codon[x] if 2-x != rf else convert for x in [0, 1, 2]))
                aa_new = new_codon.translate()

                mtype = self.get_AAconversion_type(codon1 = ref_codon, codon2 = new_codon)

                ### If conversion leads to Ref change or REf change keep
                codonpos = codon_start
                ctype = mtype
                ref = f"{ref_codon}>{aa_ref}"
                alt = f"{alt_codon}>{aa_alt}"
                convert = f"{new_codon}>{aa_new}"
                self.add_BEguides(name, guide, pam, strand, snvpos, score, start, end,ref,codonpos,alt,convert,ctype,bystander,n)
                n += 1


    def extract_Seqs(self, SNV_pos):
        genes = SeqIO.index(self.fasta_path, "fasta")
        c = 'M' if self.chrom == 'MT' else self.chrom
        self.extracted_seq = str(genes[f"chr{c}"][SNV_pos- self.search_window:SNV_pos + self.search_window].seq)
        #self.extracted_seq = str(genes.seq["chr{}".format(c)][SNV_pos - self.search_window:SNV_pos + self.search_window])
        # replace with ref allele with variant
        self.extracted_seq = Seq(self.extracted_seq[0:self.search_window] + self.NC_alt_allele + self.extracted_seq[self.search_window + 1:]).upper()
        return self.extracted_seq

    def add_guides(self,name,guide,pam_found,strand,snvpos,score,start,end,n):
        self.guides_found['HGVS_ID'].append(self.hgvs_id),
        self.guides_found['Editor'].append(name)
        self.guides_found['Guide_ID'].append(f'{name}_{self.chrom}.{self.alleleid}.{n}{strand}')
        self.guides_found['Chr'].append(self.chrom)
        self.guides_found['Score'].append(score)
        self.guides_found['Start'].append(start)
        self.guides_found['End'].append(end)
        self.guides_found['Strand'].append(strand)
        self.guides_found['Pam'].append(pam_found)
        self.guides_found['gRNA'].append(str(guide))
        self.guides_found['SNV Position'].append(snvpos)
        self.guides_found['Variant_Molecular_Consequence'].append(self.mc)


    def add_BEguides(self,name,guide,pam_found,strand,snvpos,score,start,end,codonpos,ref,alt,convert,ctype,bystander,n):
        self.BEguides_found['HGVS_ID'].append(self.hgvs_id),
        self.BEguides_found['Guide_ID'].append(f'{name}_{self.chrom}.{self.alleleid}.{n}{strand}')
        self.BEguides_found['Base Editor'].append(name)
        self.BEguides_found['Chr'].append(self.chrom)
        self.BEguides_found['Start'].append(start)
        self.BEguides_found['End'].append(end)
        self.BEguides_found['gRNA'].append(str(guide))
        self.BEguides_found['Pam'].append(pam_found)
        self.BEguides_found['SNV Position'].append(snvpos)
        self.BEguides_found['Codon Position'].append(codonpos)
        self.BEguides_found['Strand'].append(strand)
        self.BEguides_found['Reference (Codon>AA)'].append(ref)
        self.BEguides_found['Alternate (Codon>AA)'].append(alt)
        self.BEguides_found['BE Converted (Codon>AA)'].append(convert)
        self.BEguides_found['Conversion Type'].append(ctype)
        self.BEguides_found['Bystander'].append(bystander)
        self.BEguides_found['Variant_Molecular_Consequence'].append(self.mc)



    def find_pamLast(self,sitelen,target_start,guidelen,strand,win_size,guide_temp):
        guide, snvpos, start, end = None,None,None,None

        # check if in targetable window
        if target_start <= (self.search_window - win_size[0]) and target_start >= (self.search_window - win_size[1]):
            guide = guide_temp[target_start:target_start + guidelen]

            if strand == '+':
                start = (self.SNV_chr_pos - self.search_window) + target_start
                end = (self.SNV_chr_pos - self.search_window) + target_start + sitelen

            else:
                end = self.SNV_chr_pos + (self.search_window - target_start)
                start = self.SNV_chr_pos + (self.search_window - target_start) - sitelen


            snvpos = (self.search_window - target_start) if strand == "-" else (self.search_window - target_start) - 1
        return guide,snvpos,start,end

    def find_pamFirst(self,i,sitelen,target_start,guidelen,strand,win_size,guide_temp):
        guide, snvpos, start, end = None,None,None,None

        # check if in targetable window
        if target_start <= (self.search_window + win_size[1]) and target_start >= (self.search_window + win_size[0]):
            guide = guide_temp[target_start - guidelen:target_start]

            if strand == '+':
                start = (self.SNV_chr_pos - self.search_window) + i
                end = start + target_start
            else:
                end = (self.SNV_chr_pos - self.search_window) + i
                start = end + target_start
            snvpos =  (target_start - self.search_window) if strand == "-" else (target_start - self.search_window) - 1

        return guide,snvpos,start,end

    def get_guide_set(self,name,pam, pamISfirst, win_size, scoring, guidelen,BEmode):
        '''
        :param pam: pam seq ex:'NGG'
        param pamISfirst: 5'or3'PAM ex:True/False
        :param win_size: list containing upper and lower limits of the targetable window. Ex:[4,8]
        :param score: Boolean Optional Deonch scoring used for spCas9 only
        :param search window: intial search + or - SNV site
        :param guidelen: guide without pam length
        :return: Guide Dictionary
        '''
        guides =[]

        pamlen = len(pam)
        sitelen = guidelen + pamlen

        if len(self.extracted_seq) == 0: #if a extracted sequence is not already set, set it
            self.extracted_seq = self.extract_Seqs(SNV_pos=self.SNV_chr_pos)
        n = 1
        #Narrow based on guide params
        for strand in ["-","+"]:
            print(f"ERRO BIZARRO: 'ValueError: Mixed RNA/DNA found' -> {self.extracted_seq}")
            search_seq = self.extracted_seq if strand == "+" else self.extracted_seq.reverse_complement()
            guide_temp = search_seq
            pam_index = SeqUtils.nt_search(str(search_seq), pam)[1:]
            print('pam index',pam_index)
            print('strand',strand)
            print(pam,search_seq)

            for i in pam_index:

                if pamISfirst == False:# 3' PAM
                    target_start = i - guidelen
                    guide,snvpos,start,end = self.find_pamLast(sitelen,target_start,guidelen,strand,win_size,guide_temp)
                    # print(guide,snvpos,start,end)

                    if guide != None:
                        pam_found = str(guide_temp[i:i + pamlen])
                        if scoring == 'doench':
                            score = self.calcDoenchScores(guide_temp[target_start -3:target_start + sitelen + 4])
                        else:
                            score = '-'
                        guides.append([name, guide, pam_found, strand, snvpos, score,start, end])

                        if BEmode == False:
                            self.add_guides(name,guide,pam_found,strand,snvpos,score,start,end,n)
                            n += 1

                else:
                    target_start = i + sitelen
                    guide, snvpos, start, end = self.find_pamFirst(i, sitelen, target_start, guidelen, strand, win_size,guide_temp)

                    if guide != None:
                        pam_found = str(guide_temp[i:i + pamlen])
                        score = '-'

                        guides.append([name, guide, pam_found, strand, snvpos, score, start, end])

                        if BEmode == False:
                            self.add_guides(name, guide, pam_found, strand, snvpos, score, start, end,n)
                            n += 1

        return guides

    def get_Guides(self, search_params, BEsearch_params = None):

        for name, params, in search_params.items():
            scoring = 'doench' if name == 'spCas9' else None
            pam, pamISfirst, guidelen, win_size = params[0:4]
            guides = self.get_guide_set(name,pam, pamISfirst, win_size, scoring, guidelen,BEmode = False)


        # if BE mode is on
        if BEsearch_params != None:

            for k, params, in BEsearch_params.items():
                scoring = None
                pam, pamISfirst, guidelen, win_size = params[0][0], params[0][1], params[0][2], params[0][3]
                bguides = self.get_guide_set(k,pam, pamISfirst, win_size, scoring, guidelen,BEmode = True)

                # if guides are found sep neg and pos strand guides
                if len(bguides) > 0:
                    pos_guides, neg_guides = [], []
                    for g in bguides:
                        if g[3] == '+':
                            pos_guides += [g]
                        else:
                            neg_guides += [g]


                    #See if SNV can be BE edited
                    for p in range(1,len(params[1:])+1):
                        conversion = params[p][0]  # 'CT'
                        name = ",".join(
                            [n for n in params[p][1:]])  # ('BE1', 'BE2', 'BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-Gam')

                        if self.NC_alt_allele == conversion[0]:
                            if len(pos_guides) > 0:
                                self.getBE(guides=pos_guides, conversion=conversion, win_size =win_size,name = name)

                        if self.NC_alt_allele == str(Seq(conversion[0]).complement()):
                            if len(neg_guides) > 0:
                                self.getBE(guides=neg_guides, conversion=conversion,win_size =win_size,name=name)

        return self.guides_found, self.BEguides_found


'''
 test
 
datadir = "/groups/clinical/projects/editability/tables/"
fasta_path ="/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa"

search_params = {'spCas9':('NGG', False,20,[4,8], 'Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
                              'saCas9-20':('NNGRRT',False,20,[4,8],'Cas9 S. Aureus 21 base guide'),
                              'CasX':('TTCN',True,20,[4,8],'Cas12e'),
                              'Cpf1':('TTTV',True,23,[4,8],'TTT(A/C/G)-23bp - Cas12a (Cpf1)')
                              }

BE_search_params = {'spCas9-def': [('NGG', False, 20, [4, 8]), ('CT','BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-Gam'),('AG','ABE7.9', 'ABE7.10', 'ABEmax')]}
qtype= 'hgvs'
queries = ['NM_000532.5(PCCB):c.1316A>G (p.Tyr439Cys)', 'NM_000518.5(HBB)c.114G>A', 'NM_000517.6(HBA2):c.99G>A', 'NM_005886.2(KATNB1)c.1A>G']
variantseq_dict = get_seqinfo(queries,qtype,datadir)  

for hgvs_id, data in variantseq_dict.items():

    dh = DataHandler(hgvs_id,data,fasta_path)
    guides_found, BEguides_found = dh.get_Guides(search_params,BE_search_params)
    for k,v in guides_found.items():
        print(k,v)
    for k,v in BEguides_found.items():
        print(k,v)


'''

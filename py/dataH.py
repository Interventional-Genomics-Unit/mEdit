import pandas as pd
import regex as re
import numpy as np
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio import SeqUtils
import math


class Variant_DataHandler:
    '''
    To call tables based on user query and query type.
    Guide Scoring, ranking and Metrics are handled in a seperate parent class
    '''

    def __init__(self):

        # paths
        self.datadir = "/groups/clinical/projects/editability/tables"
        self.processed_tables = f"{self.datadir}/processed_tables/"  # folder with clinvar tabs
        self.HGVSlookup_path = f"{self.processed_tables}HGVSlookup.csv"  # Chrom to refID key tab

        # tables and dictionaries
        self.vardf = pd.DataFrame() #output clinical df

        # Variables
        self.chrom = str()
        self.hgvs_id = str()
        self.prefix = str() # prefix of hgvs_id/query term
        self.SNV_chr_pos = int() #chrom position of variant
        self.NC_ref_allele = str()  # genomic ref allele
        self.NC_alt_allele = str()  # genomic alt allele
        self.NM_ref_allele = str()  # coding_seq ref allele
        self.NM_alt_allele = str()  # coding_seq alt allele

    def get_HGVStable(self):
        '''
        Get HGVS Lookup table - a table to find the chrom a HGVS ID resides
        '''
        hgvs_tab = pd.read_csv(self.HGVSlookup_path)
        return hgvs_tab

    def get_Chrom(self, hgvs_id):
        '''
        Search HGVS Table for Chromosome to determine which Clinvar file is needed
        '''

        self.hgvs_id = hgvs_id
        hgvs_tab = self.get_HGVStable()
        hgvs_list = hgvs_tab["HGVS"]

        if hgvs_id.startswith('NM'):
            try:
                prefix = re.search("(NM_\d+).", hgvs_id).groups()[0]
                chrom = hgvs_tab['Chr'].iloc[list(hgvs_list).index(prefix)]
                self.chrom, self.prefix = chrom, prefix
            except:
                print("No HGVSs were found for this search")

        return self.prefix, self.chrom

    def get_ClinVartable(self, query, qtype = 'HGVS'):

        if qtype == 'HGVS':
            hgvs_id = query
            self.get_Chrom(hgvs_id)
            if len(self.chrom) > 0:
                df = pd.read_csv(f"{self.processed_tables}{self.chrom}_variant.txt", dtype=object)
                df = df.loc[df["HGVS_ID"].str.startswith(self.prefix)]
                try:
                    var_desc = re.search("[0-9c.]?[.+-]*\d+\w+>\w+", self.hgvs_id).captures()[0]
                    df = df.loc[df["HGVS_ID"].str.contains(var_desc)]
                    self.vardf = df
                    self.NC_ref_allele, self.NC_alt_allele = df.RefAlleleVCF.iloc[0], df.AltAlleleVCF.iloc[
                        0]

                    self.SNV_chr_pos = int(df.PositionVCF.iloc[0])
                    if df['Coding Strand'].iloc[0] == '-':
                        self.NM_ref_allele, self.NM_alt_allele = str(Seq(self.NC_ref_allele).complement()), str(Seq(self.NC_alt_allele).complement())
                    else:
                        self.NM_ref_allele, self.NM_alt_allele =self.NC_ref_allele, self.NC_alt_allele


                except IndexError:
                    print(f"No HGVSs were found {hgvs_id} this search")
        if qtype == 'Coord':
            # TODO: implment search by coordinats setting
            pass

        return self.vardf



class DataHandler(Variant_DataHandler):
    '''
    search hg38 fasta for guides given reqs
    pam, len, targe table window and optionally scores
    '''

    def __init__(self):
        super().__init__()

        #paths
        self.hg38_path ="/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa"

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
        self.guides_found ={'HGVS_ID':[],'Name': [],'Chr': [],
                            'Start': [],'End': [],
                            'Strand': [],'gRNA': [],
                            'Pam': [],'Score':[],
                            'SNV Position': []}

        self.BEguides = {'HGVS_ID':[],
                         'Name':[],
                         'Chr':[],
                         'Start':[],
                         'End':[],
                         'gRNA': [],
                         'Pam' : [],
                         'SNV Position': [],
                         'Codon Position':[],
                         'Strand':[],
                         'Reference (Codon>AA)':[],
                         'Alternate (Codon>AA)': [],
                         'BE Converted (Codon>AA)': [],
                         'Conversion Type':[],
                         'Bystander': []
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
    def get_AAconversion_type(self,codon1,codon2):
        '''
        codon1: codon of Alt allele to be changed by BE
        codon2: codon after changed by BE
        '''
        aa_groups = [["G","A","V","L","I"],
                     ["S","C","U","T","M"],
                     ["F","Y","W"],
                     ["H","K","R"],
                     ["D","E","N","Q"]]
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

    def getBE(self,guides,conversion,win_size):
        '''
        Finds codon level SNV and determines if the Base Editor Conversion can work
        *TO DO : Check HGVS term to find protein change and then bypass this 'manual' conversion

        '''

        BE_guides = {'HGVS_ID':[],
                     'Name':[],
                     'Chr':[],
                     'Start':[],
                     'End':[],
                     'gRNA': [],
                     'Pam' : [],
                     'SNV Position': [],
                     'Codon Position':[],
                     'Strand':[],
                     'Reference (Codon>AA)':[],
                     'Alternate (Codon>AA)': [],
                     'BE Converted (Codon>AA)': [],
                     'Conversion Type':[],
                     'Bystander': []
                     }

        f_delta = conversion # CT or AG
        r_delta = [str(Seq(conversion[0]).complement()),str(Seq(conversion[1]).complement())] ## the complment so GA or TC

        #Determine Reading frame
        coding_type, exon_pos = self.codingType()
        coding_strand = self.vardf['Coding Strand'].iloc[0]
        guide_strand = guides['Strand'][0]
        ref_allele = self.NM_ref_allele
        alt_allele = self.NM_alt_allele
        mtype = ''

        for i in range(len(guides['gRNA'])):

            gseq = Seq(guides['gRNA'][i])
            SNVpos = guides['SNV Position'][i]
            target_bases = Seq(gseq[win_size[0] - 1:win_size[1] + 1])  # Bases inside 4-8 window

            # Converted case
            convert = str(f_delta[1] if alt_allele == f_delta[0] else r_delta[1])
            bystander = target_bases.count(f_delta[0]) - 1

            if coding_type != 'exon':
                alt = alt_allele
                ref = ref_allele

            else:
                ## In Exon
                #Determine codon and tranlated product
                coding_pos = round(int(exon_pos) /3,1) #find exon positiom
                if guide_strand == coding_strand:
                    rf = 2 if coding_pos - round(coding_pos) == 0 else 1 if (coding_pos - round(
                        coding_pos)) == 0.5 else 0
                else:
                    rf = 0 if coding_pos - round(coding_pos) == 0 else 1 if (coding_pos - round(
                        coding_pos)) == 0.5 else 2

                codon_start = int(int(SNVpos) - rf)

                #Alternative codon
                alt_codon = Seq(gseq[codon_start:codon_start+3])
                if coding_strand != guide_strand:
                    alt_codon = alt_codon.reverse_complement()
                aa_alt = alt_codon.translate()

                #Reference codon
                ref_codon = Seq("".join(alt_codon[x] if x != rf else ref_allele for x in [0,1,2]))
                if coding_strand != guide_strand:
                    ref_codon = Seq("".join(alt_codon[x] if 2-x != rf else ref_allele for x in [0, 1, 2]))
                aa_ref = ref_codon.translate()

                ##Converted allele
                convert = f_delta[1] if alt_allele == f_delta[0] else r_delta[1]
                new_codon = Seq("".join(alt_codon[x] if x != rf else convert for x in [0,1,2]))
                if coding_strand != guide_strand:
                    new_codon = Seq("".join(alt_codon[x] if 2-x != rf else convert for x in [0, 1, 2]))
                aa_new = new_codon.translate()

                mtype = self.get_AAconversion_type(codon1 = ref_codon, codon2 = new_codon)

                ### If conversion leads to Ref change or REf change keep

            if coding_type != 'exon' or mtype == 'Synonymous' or mtype == 'Silent' or mtype == 'Conservative':
                BE_guides['HGVS_ID'] += [self.hgvs_id],
                BE_guides['Name'] += guides['Name'][i],
                BE_guides['Chr'] += guides['Chr'][i],
                BE_guides['Start'] += guides['Start'][i],
                BE_guides['End'] += guides['End'][i],
                BE_guides['gRNA'] += [str(gseq)]
                BE_guides['Pam']+= [guides['Pam'][i]]
                BE_guides['SNV Position']+= [SNVpos]
                BE_guides['Strand'] += [guide_strand]

                if coding_type != 'exon':
                    BE_guides['Codon Position'] += ['NA']
                    BE_guides['Reference (Codon>AA)'] += [ref]
                    BE_guides['Alternate (Codon>AA)'] += [alt]
                    BE_guides['BE Converted (Codon>AA)'] += [convert]
                    BE_guides['Conversion Type'] += ['NA']
                    BE_guides['Bystander'] += [bystander]
                else:
                    BE_guides['Codon Position'] += [codon_start]
                    BE_guides['Reference (Codon>AA)'] += [f"{ref_codon}>{aa_ref}"]
                    BE_guides['Alternate (Codon>AA)'] += [f"{alt_codon}>{aa_alt}"]
                    BE_guides['BE Converted (Codon>AA)'] += [f"{new_codon}>{aa_new}"]
                    BE_guides['Conversion Type'] += [mtype]
                    BE_guides['Bystander'] += [bystander]

        return BE_guides

    def extract_Seqs(self, SNV_pos):
        genes = Fasta(self.hg38_path)
        self.extracted_seq = str(genes[f"chr{self.chrom}"][SNV_pos- self.search_window:SNV_pos + self.search_window])
        # replace with ref allele with variant
        self.extracted_seq = Seq(self.extracted_seq[0:self.search_window] + self.NC_alt_allele + self.extracted_seq[self.search_window + 1:])
        return self.extracted_seq

    def add_guides(self,guides):

        for k,v in guides.items():
            self.guides_found[k] += v

        return self.guides_found
    def add_BEguides(self,beguides):

        for k,v in beguides.items():
            self.BEguides[k] += v

        return self.BEguides

    def get_Guides(self, search_params, BEsearch_params = None):
        '''
                search_params = {'spCas9':('NGG', False,20,[4,8], 'Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
                              'saCas9-20':('NNGRRT',False,20,[4,8],'Cas9 S. Aureus 21 base guide'),
                              'CasX':('TTCN',True,20,[4,8],'Cas12e'),
                              'Cpf1':('TTTV',True,23,[4,8],'TTT(A/C/G)-23bp - Cas12a (Cpf1)')
                              }

                BEsearch_params = {'spCas9-def': [('NGG', False, 20, [4, 8]),
                ('CT', 'BE3', 'BE4', 'BE4max', 'BE4-Gam'), ('AG', 'ABE7.9', 'ABE7.10', 'ABEmax')]}

        '''

        for name, params, in search_params.items():
            scoring = 'doench' if name == 'spCas9' else None
            pam, pamISfirst, guidelen, win_size = params[0], params[1], params[2], params[3]

            guides = self.get_guide_set(name,pam, pamISfirst, win_size, scoring, guidelen)

            if len(guides['gRNA']) > 0:

                self.add_guides(guides)


        # if BE mode is on
        if BEsearch_params != None:

            ct = 1
            for k, params, in BEsearch_params.items():
                name = ",".join(
                    [n for n in params[ct][1:]])  # ('BE3', 'HF-BE3', 'BE4', 'BE4max')
                scoring = None
                pam, pamISfirst, guidelen, win_size = params[0][0], params[0][1], params[0][2], params[0][3]
                bguides = self.get_guide_set(name,pam, pamISfirst, win_size, scoring, guidelen)
                ct+=1
                # if guides are found sep neg and pos strand guides
                if len(bguides['gRNA']) > 0:

                    mat = np.matrix(list(bguides.values()))
                    pos_guides, neg_guides = {}, {}

                    # negative strand guides
                    try:
                        sel = np.where(mat[list(bguides.keys()).index('Strand')] == '-')[1]
                        k = list(bguides.keys())
                        for i in range(len(k)):
                            neg_guides[k[i]] = mat[i, sel].tolist()[0]
                    except:
                        pass

                    # positive strand guides
                    try:
                        sel = np.where(mat[list(bguides.keys()).index('Strand')] == '+')[1]
                        k = list(bguides.keys())
                        for i in range(len(k)):
                            pos_guides[k[i]] = mat[i, sel].tolist()[0]
                    except:
                        pass

                    #See if SNV can be BE edited
                    for p in range(1,len(params[1:])+1):

                        try:
                            conversion = params[p][0]  # 'CT'
                            name = ",".join(
                                [n for n in params[p][1:]])  # ('BE1', 'BE2', 'BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-Gam')

                            if self.NC_alt_allele == conversion[0]:
                                if len(pos_guides.keys()) > 0:
                                    beguides = self.getBE(guides=pos_guides, conversion=conversion, win_size =win_size )

                            if self.NC_alt_allele == str(Seq(conversion[0]).complement()):
                                if len(neg_guides.keys()) > 0:
                                    beguides = self.getBE(neg_guides, conversion=conversion,win_size =win_size)


                            if len(beguides['gRNA']) > 0:
                                self.add_BEguides(beguides)
                        except:
                            pass


            return self.guides_found, self.BEguides
        else:
            return self.guides_found


    def get_guide_set(self,name,pam, pamISfirst, win_size, scoring, guidelen):
        '''
        :param pam: pam seq ex:'NGG'
        param pamISfirst: 5'or3'PAM ex:True/False
        :param win_size: list containing upper and lower limits of the targetable window. Ex:[4,8]
        :param score: Boolean Optional Deonch scoring used for spCas9 only
        :param search window: intial search + or - SNV site
        :param guidelen: guide without pam length
        :return: Guide Dictionary
        '''
        guides = {
            'HGVS_ID': [],
            'Name': [],
            'Chr': [],
            'Start': [],
            'End': [],
            'Strand': [],
            'gRNA': [],
            'Pam': [],
            'Score':[],
            'SNV Position': []}

        pamlen = len(pam)
        sitelen = guidelen + pamlen

        if len(self.extracted_seq) == 0: #if a extracted sequence is not already set, set it
            self.extracted_seq = self.extract_Seqs(SNV_pos=self.SNV_chr_pos)


        #Narrow based on guide params
        cnt = 0
        for strand in ["-","+"]:
            search_seq = self.extracted_seq if strand == "+" else self.extracted_seq.reverse_complement()
            guide_temp = search_seq

            pam_index = SeqUtils.nt_search(str(search_seq), pam)[1:]

            for i in pam_index:
                if pamISfirst == False:

                    target_start = i - guidelen

                    if target_start <= (self.search_window - win_size[0]) and target_start >= (self.search_window - win_size[1]):

                        guide = guide_temp[target_start:target_start + guidelen]

                        if scoring == 'doench':
                            guides['Score'].append(self.calcDoenchScores(guide_temp[target_start -3:target_start + sitelen + 4]))
                        else:
                            guides['Score'].append('-')


                        if strand == '+':
                            start = (self.SNV_chr_pos - self.search_window) + target_start
                            end = (self.SNV_chr_pos - self.search_window) + target_start + sitelen

                        else:
                            end = self.SNV_chr_pos + (self.search_window - target_start)
                            start = self.SNV_chr_pos + (self.search_window - target_start) - sitelen
                        guides['HGVS_ID'].append(self.hgvs_id),
                        guides['Name'].append(f'{name}_{self.chrom}{start}{strand}')
                        guides['Chr'].append(self.chrom)
                        guides['Start'].append(start)
                        guides['End'].append(end)
                        guides['Strand'].append(strand)
                        guides['Pam'].append(str(guide_temp[i:i+pamlen]))
                        guides['gRNA'].append(str(guide))
                        guides['SNV Position'].append(
                            (self.search_window - target_start) if strand == "-" else (self.search_window - target_start) - 1)
                else:
                    target_start = i + sitelen

                    if target_start <= (self.search_window + win_size[1]) and target_start >= (self.search_window + win_size[0]):

                        guide = guide_temp[target_start - sitelen:target_start]

                        guides['Score'].append('-')

                        if strand == '+':
                            start = (self.SNV_chr_pos - self.search_window) + i
                            end = start + target_start
                        else:
                            end = (self.SNV_chr_pos - self.search_window) + i
                            start = end + target_start
                        guides['HGVS_ID'].append(self.hgvs_id),
                        guides['Name'].append(f'{name}_{self.chrom}{start}{strand}')
                        guides['Chr'].append(self.chrom)
                        guides['End'].append(end)
                        guides['Start'].append(end)
                        guides['Strand'].append(strand)
                        guides['Pam'].append(str(guide_temp[i:i + pamlen]))
                        guides['gRNA'].append(str(guide))
                        guides['SNV Position'].append(
                            (target_start - self.search_window) if strand == "-" else (target_start - self.search_window) - 1)
        return guides


'''
 test
search_params = {'spCas9':('NGG', False,20,[4,8], 'Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
                              'saCas9-20':('NNGRRT',False,20,[4,8],'Cas9 S. Aureus 21 base guide'),
                              'CasX':('TTCN',True,20,[4,8],'Cas12e'),
                              'Cpf1':('TTTV',True,23,[4,8],'TTT(A/C/G)-23bp - Cas12a (Cpf1)')
                              }

BEsearch_params = {'spCas9-def': [('NGG', False, 20, [4, 8]), ('CT', 'BE1', 'BE2', 'BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-Gam'),
 ('AG','ABE7.9', 'ABE7.10', 'ABEmax')]}
hgvs_id = 'NM_000152.5(GAA):c.271G>T'
qtype= 'HGVS'

dh = DataHandler()
clindf = dh.get_ClinVartable('NM_000016.6(ACADM):c.727C>T (p.Arg243Ter)',qtype)
guides_found, BEguides_found = dh.get_Guides(fg.search_params,fg.BE_search_params)
for k,v in guides_found.items():
    print(k,v)
for k,v in BEguides_found.items():
    print(k,v)

---output----

Chr ['17', '17', '17', '17', '17', '17']
Start [80104842, 80104841, 80104840, 80104839, 80104838, 80104852]
End [80104865, 80104864, 80104863, 80104862, 80104861, 80104875]
Strand ['-', '-', '-', '-', '-', '+']
gRNA ['GGCGCAAACGAAGCGGCTGT', 'GCGCAAACGAAGCGGCTGTT', 'CGCAAACGAAGCGGCTGTTG', 'GCAAACGAAGCGGCTGTTGG', 'CAAACGAAGCGGCTGTTGGG', 'CTTCGTTTGCGCCCCTGACA']
Pam ['TGG', 'GGG', 'GGG', 'GGG', 'GGG', 'AGG']
Score [6, 0, 7, 13, 7, 36]
SNV Position [8, 7, 6, 5, 4, 4]


{'ABE7.9,ABE7.10,ABEmax': {'Chr': ['17', '17', '17', '17', '17'],
 'Start': ['80104842', '80104841', '80104840', '80104839', '80104838'], 
 'End': ['80104865', '80104864', '80104863', '80104862', '80104861'], 
 'gRNA': ['GGCGCAAACGAAGCGGCTGT', 'GCGCAAACGAAGCGGCTGTT', 'CGCAAACGAAGCGGCTGTTG', 'GCAAACGAAGCGGCTGTTGG', 'CAAACGAAGCGGCTGTTGGG'],
  'Pam': ['TGG', 'GGG', 'GGG', 'GGG', 'GGG'],
   'SNV Position': ['8', '7', '6', '5', '4'],
    'Codon Position': [6, 5, 4, 3, 2], 
    'Strand': ['-', '-', '-', '-', '-'], 
    'Reference (Codon>AA)': ['GTT>V', 'GTT>V', 'GTT>V', 'GTT>V', 'GTT>V'],
     'Alternate (Codon>AA)': ['GTT>V', 'GTT>V', 'GTT>V', 'GTT>V', 'GTT>V'], 
     'BE Converted (Codon>AA)': ['CTT>L', 'CTT>L', 'CTT>L', 'CTT>L', 'CTT>L'], 
     'Conversion Type': ['Conservative', 'Conservative', 'Conservative', 'Conservative', 'Conservative'],
      'Bystander': [2, 2, 3, 3, 2]}}

'''

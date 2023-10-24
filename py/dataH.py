from Bio.Seq import Seq
from Bio import SeqUtils
import math


class DataHandler:
    """
    search for guides given a genomic sequence and SNV info
    """

    def __init__(self, query, strand, ref, alt, feature_annotation, extracted_seq, rf, coord):
        """
        :param query: ex: 'NM_000532.5(PCCB):c.1316A>G (p.Tyr439Cys)' or 'chr19:136327650A>G'
        :param strand: ex. '-' or '+
        :param ref: ex. 'A'
        :param alt: ex. 'G'
        :param feature_annotation: ex: 'exon'
        :param extracted_seq: 'CCCACAGGGCCCTCACCTGCAGATTGTGATTGTGGCCGCACAGGTAGGCAGTGACCCCGT'
        :param rf : ex: '2'
        :param coord : ex: 'chr19:136327650'
        """
        # search data
        self.NC_ref_allele = str(ref).upper()
        self.NC_alt_allele = str(alt).upper()
        self.strand = strand  # coding_strand
        self.NM_ref_allele = self.NC_ref_allele if self.strand == '+' else str(Seq(self.NC_ref_allele).complement())
        self.NM_alt_allele = self.NC_alt_allele if self.strand == '+' else str(Seq(self.NC_alt_allele).complement())
        self.SNV_chr_pos = int(coord.split(':')[1])
        self.query = query
        self.rf = rf
        self.extracted_seq = str(extracted_seq)
        self.annotation = feature_annotation
        self.coord = coord
        self.chrom = coord.split(':')[0].replace('chr', '')

        # search params
        self.pam = str()  # Ex. 'NGG'
        self.pamISfirst = False  # Boolean
        self.win_size = list()  # Ex. list [4,8]
        self.scoring = None  # Ex. True/False
        self.guidelen = 20

        # outputs
        self.guides_found = {'QueryTerm': [], 'Editor': [], 'Guide_ID': [], 'Coordinates': [],
                             'Strand': [], 'gRNA': [],
                             'Pam': [], 'Doench Score': [],
                             'SNV Position': [], 'Ref>Alt': [], 'Annotation': []}

        self.BEguides_found = {'QueryTerm': [], 'Base Editor': [], 'Guide_ID': [], 'Coordinates': [],
                               'gRNA': [], 'Pam': [], 'SNV Position': [], 'Strand': [],
                               'Reference (Codon>AA)': [], 'Alternate (Codon>AA)': [], 'BE Converted (Codon>AA)': [],
                               'Conversion Type': [], 'Bystander': [], 'Annotation': []}

        # tables and dictionaries
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

    def set_guide_search_params(self, pam, pamISfirst, win_size, scoring, guidelen):
        self.pam = pam
        self.pamISfirst = pamISfirst
        self.win_size = win_size
        self.scoring = scoring
        self.guidelen = guidelen

    def calcDoenchScores(self, seq):
        """
        Input is a 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'
        """
        global gcWeight
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

    def find_codon(self, snv_rel_pos):
        codon = self.extracted_seq[int(snv_rel_pos - self.rf): int((snv_rel_pos - self.rf) + 3)]
        if self.strand == '-':
            codon = Seq(codon).complement()
        return codon

    @staticmethod
    def get_AAconversion_type(codon1, codon2):
        """
        codon1: codon of Alt allele to be changed by BE
        codon2: codon after changed by BE
        """
        aa_groups = [["G", "A", "V", "L", "I", "M", "F", "Y", "W"],
                     ["S", "Q", "T", "N"],
                     ["C", "G", "P"],
                     ["D", "E"],
                     ["K", "H", "R", "Q"]]
        codon1, codon2 = Seq(codon1), Seq(codon2)
        aa1, aa2 = codon1.translate(), codon2.translate()
        mtype = ""
        if aa1 == aa2:
            mtype = 'Synonymous'
            if codon1 == codon2:
                mtype = 'Silent'
        else:
            if codon2 in ['TAA', 'TAG', 'TGA']:
                mtype = 'Nonsense'
            elif ([aa_groups.index(x) for x in aa_groups if str(aa1) in x] ==
                  [aa_groups.index(x) for x in aa_groups if str(aa2) in x]):
                mtype = 'Conservative'
            else:
                mtype = 'Non-conservative'
        return mtype

    def getBE(self, guides, conversion, win_size, name):
        """
        Finds codon level SNV and determines if the Base Editor Conversion can work
        """
        coding_strand = self.strand
        snv_rel_pos = len(self.extracted_seq) / 2

        for i in range(len(guides)):
            editor, guide, pam_found, guide_strand, snvpos, score, start, end = guides[i]

            target_bases = Seq(guide[win_size[0] - 1:win_size[1] + 1])  # Bases inside 4-8 window

            # Converted case
            convert = str(conversion[1])
            bystander = target_bases.count(conversion[0]) - 1

            if self.annotation != 'exon':
                ctype = 'NA'
                self.add_BEguides(name,
                                  guide,
                                  pam_found,
                                  guide_strand,
                                  snvpos,
                                  start,
                                  end,
                                  self.NM_ref_allele,
                                  self.NM_alt_allele,
                                  convert,
                                  ctype,
                                  bystander)

            else:
                ## In Exon
                #Determine codon and translated product
                alt_codon = self.find_codon(snv_rel_pos)

                # Alternative codon
                if coding_strand != guide_strand:
                    alt_codon = Seq(alt_codon).reverse_complement()
                print("\n===============\nThe method 'translate' might need an input that is currently not provided\n===========\n")
                aa_alt = alt_codon.translate()

                # Reference codon
                ref_codon = Seq("".join(alt_codon[x] if x != abs(self.rf) else self.NM_ref_allele for x in [0,1,2]))
                aa_ref = ref_codon.translate()

                ##Converted allele
                convert = convert if self.NM_alt_allele == conversion[0] else str(Seq(convert).complement())

                new_codon = Seq("".join(alt_codon[x] if x != abs(self.rf) else convert for x in [0,1,2]))
                aa_new = new_codon.translate()

                mtype = self.get_AAconversion_type(codon1=ref_codon, codon2=new_codon)

                ### If conversion leads to Ref change or REf change keep
                ctype = mtype
                ref = f"{ref_codon}>{aa_ref}"
                alt = f"{alt_codon}>{aa_alt}"
                convert = f"{new_codon}>{aa_new}"
                self.add_BEguides(name,
                                  guide,
                                  pam_found,
                                  guide_strand,
                                  snvpos,
                                  start,
                                  end,
                                  ref,
                                  alt,
                                  convert,
                                  ctype,
                                  bystander)

    def add_guides(self, name, guide, pam_found, strand, snvpos, score, start, end):
        self.guides_found['QueryTerm'].append(self.query)
        self.guides_found['Editor'].append(name)
        self.guides_found['Guide_ID'].append(f'{name}_')
        self.guides_found['Coordinates'].append(f'{self.chrom}:{start}-{end}')
        self.guides_found['Doench Score'].append(score)
        self.guides_found['Strand'].append(strand)
        self.guides_found['Pam'].append(str(pam_found))
        self.guides_found['gRNA'].append(str(guide))
        self.guides_found['Ref>Alt'].append(f"{self.NM_ref_allele}>{self.NM_alt_allele}")
        self.guides_found['SNV Position'].append(snvpos)
        self.guides_found['Annotation'].append(self.annotation)

    def add_BEguides(self, name, guide, pam_found, strand, snvpos, start, end, ref, alt, convert, ctype, bystander):
        self.BEguides_found['QueryTerm'].append(self.query)
        self.BEguides_found['Guide_ID'].append(f'{name}_')
        self.BEguides_found['Base Editor'].append(name)
        self.BEguides_found['Coordinates'].append(f'{self.chrom}:{start}-{end}')
        self.BEguides_found['gRNA'].append(str(guide))
        self.BEguides_found['Pam'].append(str(pam_found))
        self.BEguides_found['SNV Position'].append(snvpos)
        self.BEguides_found['Strand'].append(strand)
        self.BEguides_found['Reference (Codon>AA)'].append(ref)
        self.BEguides_found['Alternate (Codon>AA)'].append(alt)
        self.BEguides_found['BE Converted (Codon>AA)'].append(convert)
        self.BEguides_found['Conversion Type'].append(ctype)
        self.BEguides_found['Bystander'].append(bystander)
        self.BEguides_found['Annotation'].append(self.annotation)

    def get_guide_set(self, name, pam, pamISfirst,scoring, win_size, guidelen, BEmode):
        """
        :param name:
        :param pam: pam seq ex:'NGG'
        :param pamISfirst: 5'or3'PAM ex:True/False
        :param win_size: list containing upper and lower limits of the targetable window. Ex:[4,8]
        :param score: Boolean Optional Deonch scoring used for spCas9 only
        :param search window: intial search + or - SNV site
        :param guidelen: guide without pam length
        :return: Guide Dictionary
        """
        guides = []
        pamlen = len(pam)
        sitelen = guidelen + pamlen
        snv_rel_pos = int(len(self.extracted_seq)/2)

        if BEmode:
            win_size = [win_size[0] - guidelen,win_size[1] -guidelen]

        pam_min, pam_max = int((snv_rel_pos - win_size[1]))- 1, int((snv_rel_pos - win_size[0])) - 1

        if pamISfirst == True:
            pam_min, pam_max = pam_min + pamlen, pam_max + pamlen

        # Narrow based on guide params
        for search_strand in ["-", "+"]:
            search_seq = Seq(self.extracted_seq) if search_strand == "+" else Seq(self.extracted_seq).reverse_complement()
            pam_index = SeqUtils.nt_search(str(search_seq), pam)[1:]

            for i in pam_index:
                if i in range(pam_min, pam_max + 1):
                    if not pamISfirst:  # 3' PAM
                        target_start = i - guidelen
                        guide = search_seq[i - guidelen:i]
                        pam_found = str(search_seq[i:i + pamlen])
                        if scoring == 'doench':
                            score = self.calcDoenchScores(search_seq[target_start -3:target_start + sitelen + 4])
                        else:
                            score = '-'
                    else:
                        target_start = i + pamlen
                        guide = search_seq[target_start: i + sitelen]
                        pam_found = search_seq[i:target_start]
                        score = '-'

                    snvpos = snv_rel_pos - target_start
                    start = self.SNV_chr_pos - snvpos
                    end = start + sitelen

                    guides.append([name, guide, pam_found, search_strand, snvpos, score, start, end])

                    if not BEmode:
                        self.add_guides(name, guide, pam_found, search_strand, snvpos, score, start, end)
        return guides

    def get_Guides(self, search_params, BEsearch_params=None):
        for name, params, in search_params.items():
            scoring = 'doench' if name == 'spCas9' else None
            pam, pamISfirst, guidelen, dsb_loc = params[0:4]
            win_size = [int(dsb_loc)-7, int(dsb_loc)+7]
            guides = self.get_guide_set(name, pam, pamISfirst, scoring, win_size, guidelen, BEmode=False)

        # if BE mode is on
        if BEsearch_params is not None:
            for k, params, in BEsearch_params.items():
                scoring = None
                pam, pamISfirst, guidelen, win_size = params[0][0], params[0][1], params[0][2], params[0][3]
                bguides = self.get_guide_set(k, pam, pamISfirst, scoring, win_size, guidelen, BEmode=True)

                # if guides are found sep neg and pos strand guides
                if len(bguides) > 0:
                    pos_guides, neg_guides = [], []
                    for g in bguides:
                        if g[3] == '+':
                            pos_guides += [g]
                        else:
                            neg_guides += [g]

                    # See if SNV can be BE edited
                    for p in range(1, len(params[1:]) + 1):
                        conversion = params[p][0]  # 'CT'
                        name = ",".join(
                            [n for n in params[p][1:]])

                        if self.NC_alt_allele == conversion[0]:
                            if len(pos_guides) > 0:
                                self.getBE(guides=pos_guides, conversion=conversion, win_size=win_size, name=name)

                        if self.NC_alt_allele == str(Seq(conversion[0]).complement()):
                            if len(neg_guides) > 0:
                                self.getBE(guides=neg_guides, conversion=conversion, win_size=win_size, name=name)

        return self.guides_found, self.BEguides_found


'''
#----------------------------Test----------------------
datadir = "/groups/clinical/projects/editability/tables/"
processed_tables = "/groups/clinical/projects/editability/tables/processed_tables/"
fasta_path ="/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa.gz"

search_params= {'spCas9': ('NGG', False, 20, -2, 'Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
                'saCas9': ('NNGRRT', False,21, -2, 'Cas9 S. Aureus 21 base guide'),
                'spG': ('NGN', False, 20, -2, '20bp-NGN - SpG'),
                'SpRY-HighE': ('NRN', False,20, -2, 'High Efficiency Pam'),
                'scCas9':('NNGT',False,20,-2,'20bp-NNGT - Cas9 S. canis - high efficiency PAM, recommended'),
                'stCas9': ('NNAGAA', False,20, -2, 'Cas9 S. Thermophilus'),
                'iSpyMacCas9': ('NAA', False,20, -2, ''),
                'CasX': ('TTCN', True, 20, 18, 'Cas12e'),
                'AsCas12a': ('TTTV', True, 23, 22, 'TTT(A/C/G)-23bp - Cas12a (Cpf1)'),
                'LbCas12a': ('TTTV', True, 23, 22, 'LbCpf1'),
                'Cas12c1': ('TG', True, 23, 22, 'C2c3')}

BE_search_params = {'spCas9-def': [('NGG', False, 20, [4, 8]), ('CT', 'BE3', 'BE4', 'BE4max', 'BE4-Gam'), ('AG', 'ABE7.9', 'ABE7.10', 'ABEmax')]}

snv_info = {'11': [['NM_000518.5:c.114G>A', '-', 'C', 'T', 'exon', Seq('ATCCCCAAAGGACTCAAAGAACCTCTGGGTTCAAGGGTAGACCACCAGCAGCCTAAGGGT'), -2, 'chr11:5226778']], 
            '3': [['NM_000532.5:c.1316A>G', '+', 'A', 'G', 'exon', Seq('TGGATCTGTTTTAGGCCTATGGAGGTGCCTGTGATGTCATGAGCTCTAAGCACCTTTGTG'), 1, 'chr3:136327650']],
            '16': [['NM_000517.6:c.99G>A', '+', 'G', 'A', 'exon', Seq('CACCCCTCACTCTGCTTCTCCCCGCAGGATATTCCTGTCCTTCCCCACCACCAAGACCTA'), 2, 'chr16:173128'], 
                   ['NM_005886.3:c.1A>G', '+', 'A', 'G', 'exon', Seq('GTGGGGCTTCAGGTGCCAGCCAGCTGAAGGGTGGCCACCCCTGTGGTCACCAAGACAGCC'), 2, 'chr16:57737244']]}

for ch, data in snv_info.items():
    for d in data:
        query, strand, ref, alt, feature_annotation, extracted_seq, codons, coord = d
        dh = DataHandler(query, strand, ref, alt, feature_annotation, extracted_seq, codons, coord)
    guides_found, BEguides_found = dh.get_Guides(search_params,BE_search_params)
    for k,v in guides_found.items():
        print(k,v)
    for k,v in BEguides_found.items():
        print(k,v)


'''

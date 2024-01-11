# Installed Modules
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio.SeqUtils import seq3
# Project Modules
import scoring


class DataHandler:
    """
    search for guides given a genomic sequence and SNV info
    """

    def __init__(self, query, strand, ref, alt, feature_annotation, extracted_seq, rf, coord,gname):
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
        self.gname = gname

        # search params
        self.pam = str()  # Ex. 'NGG'
        self.pamISfirst = False  # Boolean
        self.win_size = list()  # Ex. list [4,8]
        self.gscoring = None  # Ex. True/False
        self.guidelen = 20

        # outputs
        self.guides_found = {'QueryTerm': [], 'GeneName':[],'Editor': [], 'Guide_ID': [], 'Coordinates': [],
                             'Strand': [], 'gRNA': [], 'Pam': [], 'SNV Position': [],
                             'On-Target Efficiency Score': [], 'OOF Score':[],'Ref>Alt': [], 'Annotation': []}

        self.BEguides_found = {'QueryTerm': [], 'GeneName':[],'Base Editor': [], 'Guide_ID': [], 'Coordinates': [],
                               'Strand': [],'gRNA': [], 'Pam': [], 'SNV Position': [],
                               'Hg38 Reference (Codon>AA)': [], 'Alternate (Codon>AA)': [], 'BE Converted (Codon>AA)': [],
                               'Conversion Type': [], 'Bystander': [], 'Annotation': []}


    def set_guide_search_params(self, pam, pamISfirst, win_size, gscoring, guidelen):
        self.pam = pam
        self.pamISfirst = pamISfirst
        self.win_size = win_size
        self.gscoring = gscoring
        self.guidelen = guidelen

    def find_codon(self, snv_rel_pos):
        if self.strand == '+':
            codon = self.extracted_seq[int(snv_rel_pos - self.rf): int((snv_rel_pos - self.rf) + 3)]
        else:
            adj_rf = 2 - self.rf
            codon = self.extracted_seq[int(snv_rel_pos - adj_rf): int(snv_rel_pos - adj_rf)+3]
            codon = Seq(codon).reverse_complement()
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
            editor, guide, pam_found, guide_strand, snvpos, on_score, oof_score, start, end = guides[i]

            target_bases = Seq(guide[win_size[0] - 1:win_size[1] + 1])  # Bases inside 4-8 window

            # Converted case
            convert = str(conversion[1])
            bystander = target_bases.count(conversion[0]) - 1
            if self.annotation not in ['exon','start_codon','stop_codon']:
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
                aa_alt = Seq(alt_codon).translate()

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
                aa_new = seq3(aa_new, custom_map={"*": "***"})
                aa_alt = seq3(aa_alt, custom_map={"*": "***"})
                aa_ref = seq3(aa_ref, custom_map={"*": "***"})
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

    def add_guides(self, name, guide, pam_found, strand, snvpos, on_score,oof_score,start,end):
        self.guides_found['QueryTerm'].append(self.query)
        self.guides_found['GeneName'].append(self.gname)
        self.guides_found['Editor'].append(name)
        self.guides_found['Guide_ID'].append(f'{name}_')
        self.guides_found['Coordinates'].append(f'{self.chrom}:{start}-{end}')
        self.guides_found['On-Target Efficiency Score'].append(on_score)
        self.guides_found['OOF Score'].append(oof_score)
        self.guides_found['Strand'].append(strand)
        self.guides_found['Pam'].append(str(pam_found))
        self.guides_found['gRNA'].append(str(guide))
        self.guides_found['Ref>Alt'].append(f"{self.NM_ref_allele}>{self.NM_alt_allele}")
        self.guides_found['SNV Position'].append(snvpos)
        self.guides_found['Annotation'].append(self.annotation)

    def add_BEguides(self, name, guide, pam_found, strand, snvpos, start, end, ref, alt, convert, ctype, bystander):
        self.BEguides_found['QueryTerm'].append(self.query)
        self.BEguides_found['GeneName'].append(self.gname)
        self.BEguides_found['Guide_ID'].append(f'{name}_')
        self.BEguides_found['Base Editor'].append(name)
        self.BEguides_found['Coordinates'].append(f'{self.chrom}:{start}-{end}')
        self.BEguides_found['gRNA'].append(str(guide))
        self.BEguides_found['Pam'].append(str(pam_found))
        self.BEguides_found['SNV Position'].append(snvpos)
        self.BEguides_found['Strand'].append(strand)
        self.BEguides_found['Hg38 Reference (Codon>AA)'].append(ref)
        self.BEguides_found['Alternate (Codon>AA)'].append(alt)
        self.BEguides_found['BE Converted (Codon>AA)'].append(convert)
        self.BEguides_found['Conversion Type'].append(ctype)
        self.BEguides_found['Bystander'].append(bystander)
        self.BEguides_found['Annotation'].append(self.annotation)

    def get_guide_set(self, name, pam, pamISfirst, win_size, guidelen, BEmode):
        """
        :param name:
        :param pam: pam seq ex:'NGG'
        :param pamISfirst: 5'or3'PAM ex:True/False
        :param win_size: list containing upper and lower limits of the targetable window. Ex:[4,8]
        :param search window: intial search + or - SNV site
        :param guidelen: guide without pam length
        :return: Guide Dictionary
        """
        guides = []
        pamlen = len(pam)
        sitelen = guidelen + pamlen
        snv_rel_pos = int(len(self.extracted_seq)/2)
        on_score, oof_score = '-','-'
        if BEmode:
            win_size = [win_size[0] - guidelen -1,win_size[1] -guidelen-1]

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
                        if not BEmode:
                            if pam == 'NGG' and guidelen == 20:
                                #Azmith only accurate for NGG pams
                                on_score = scoring.azimuth(cas9_sites = [str(search_seq[target_start - 3:target_start + sitelen + 4])])[0]
                            #oof_score use only for DSB
                            mh_score, oof_score = scoring.oofscore(str(search_seq[target_start - 20:target_start + sitelen + 20]))

                        pam_found = str(search_seq[i:i + pamlen])

                    else:
                        target_start = i + pamlen
                        guide = search_seq[target_start: i + sitelen]
                        pam_found = search_seq[i:target_start]
                        if 'Cas12a' in name:
                            on_score = round(scoring.deepcpf1([str(search_seq[i - 5:i + sitelen + 4])])[0][0], 2)

                    snvpos = snv_rel_pos - target_start
                    start = self.SNV_chr_pos - snvpos
                    end = start + sitelen

                    guides.append([name, guide, pam_found, search_strand, snvpos, on_score, oof_score, start, end])

                    if not BEmode:
                        self.add_guides(name, guide, pam_found, search_strand, snvpos, on_score, oof_score, start, end)
        return guides

    def get_Guides(self, search_params, BEsearch_params=None):

        for name, params, in search_params.items():
            pam, pamISfirst, guidelen, dsb_loc = params[0:4]
            win_size = [int(dsb_loc)-7, int(dsb_loc)+7]
            guides = self.get_guide_set(name, pam, pamISfirst,win_size, guidelen, BEmode=False)

        # if BE mode is on
        if BEsearch_params is not None:
            for k, params, in BEsearch_params.items():
                pam, pamISfirst, guidelen, win_size = params[0][0], params[0][1], params[0][2], params[0][3]
                bguides = self.get_guide_set(k, pam, pamISfirst, win_size, guidelen, BEmode=True)

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
                'CasX': ('TTCN', True, 20, 18, 'Cas12e'),
                'AsCas12a': ('TTTV', True, 23, 22, 'TTT(A/C/G)-23bp - Cas12a (Cpf1)')}
BE_search_params = {'spCas9-def': [('NGG', False, 20, [4, 8]), ('CT', 'BE3', 'BE4', 'BE4max', 'BE4-Gam'), ('AG', 'ABE7.9', 'ABE7.10', 'ABEmax')]}

snv_info = {'11': [['NM_000518.5:c.114G>A', '-', 'C', 'T', 'exon', Seq('ATCCCCAAAGGACTCAAAGAACCTCTGGGTTCAAGGGTAGACCACCAGCAGCCTAAGGGT'), -2, 'chr11:5226778']], 
            '3': [['NM_000532.5:c.1316A>G', '+', 'A', 'G', 'exon', Seq('TGGATCTGTTTTAGGCCTATGGAGGTGCCTGTGATGTCATGAGCTCTAAGCACCTTTGTG'), 1, 'chr3:136327650']],
            '16': [['NM_000517.6:c.99G>A', '+', 'G', 'A', 'exon', Seq('CACCCCTCACTCTGCTTCTCCCCGCAGGATATTCCTGTCCTTCCCCACCACCAAGACCTA'), 2, 'chr16:173128'], 
                   ['NM_005886.3:c.1A>G', '+', 'A', 'G', 'exon', Seq('GTGGGGCTTCAGGTGCCAGCCAGCTGAAGGGTGGCCACCCCTGTGGTCACCAAGACAGCC'), 2, 'chr16:57737244']]}

for ch, data in fg.snv_info.items():
    for d in data:
        query, tid, eid, strand, ref, alt, feature_annotation, extracted_seq, codons, coord = d
        print(f'----------{query}--------------')
        dh = DataHandler(query, strand, ref, alt, feature_annotation, extracted_seq, codons, coord)
    guides_found, BEguides_found = dh.get_Guides(fg.search_params,fg.BE_search_params)
        print(len(guides_found['gRNA'])
            #for k,v in guides_found:
    #    print(k,v)
    #for k,v in BEguides_found.items():
    #    print(k,v)

'''

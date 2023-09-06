import pandas as pd
import regex as re
from pyfaidx import Fasta
from Bio.Seq import Seq


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
        self.NC_ref_allele = str()  # genomic ref seq
        self.NC_alt_allele = str()  # genomic ref seq

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

    def get_ClinVartable(self, query, qtype):

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
                    self.SNV_chr_pos = int(df.PositionVCF)

                except IndexError:
                    print(f"No HGVSs were found {hgvs_id} this search")
        if qtype == 'Coord':
            # TODO: implment search by coordinats setting
            pass

        return self.vardf



class DataHandler(Variant_DataHandler):
    '''
    search hg38 fasta for guides given reqs
    pam, len, targetable window and optionally scores
    '''

    def __init__(self):
        super().__init__()

        #paths
        self.hg38_path ="/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa"

        #outputs
        self.guides = {}

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


    def calcDoenchScores(self, seq):
        """
        Code reproduced following paper's methods section. Thanks to Daniel McPherson for fixing it.
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

    def get_Guide(self, pam, win_size,score,guidelen=20):
        '''

        :param pam: pam seq ex:'NGG'
        :param win_size: list containing upper and lower limits of the targetable window. Ex:[4,8]
        :param score: Boolean Optional Deonch scoring used for spCas9 only
        :param guidelen: guide without pam length
        :return: Guide Dictionary
        '''
        # pos, chrom = 80104857,17
        # win_size,pam,score,guidelen = [4,8], 'NGG',True, 20
        # hgvs_id = 'NM_000152.5(GAA):c.271G>T'

        guides = {
            'Chr': [],
            'Start': [],
            'End': [],
            'Strand': [],
            'gRNA': [],
            'Pam': [],
            'Score':[],
            'SNV Position': []
        }

        genes = Fasta(self.hg38_path)
        pam_len = len(pam)
        site_len = guidelen + pam_len
        search_window = 30  # intial search + or - SNV site

        seq = str(genes[f"chr{self.chrom}"][self.SNV_chr_pos- search_window:self.SNV_chr_pos + search_window])
        # replace with variant
        seq = seq[0:search_window] + self.NC_alt_allele + seq[search_window+1 :]

        seq = Seq(seq)

        for strand in ["-","+"]:
            found, slider = 0, 0
            search_seq = seq if strand == "+" else seq.reverse_complement()
            guide_temp = seq if strand == "+" else seq.reverse_complement()

            while found != -1:
                pam_position = found + slider
                pam_start = pam_position - 1
                target_start = pam_start - guidelen

                if target_start <= (search_window - win_size[0]) and target_start >= (search_window - win_size[1]):

                    guide = guide_temp[target_start:target_start + site_len]
                    if score == True:
                        print(guide_temp[target_start - 3:target_start + site_len + 4])
                        guides['Score'].append(self.calcDoenchScores(guide_temp[target_start - 3:target_start + site_len + 4]))
                    else:
                        guides['Score'].append('-')

                    if strand == '+':
                        guides['Start'].append((self.SNV_chr_pos - search_window) + target_start)
                        guides['End'].append((self.SNV_chr_pos - search_window) + target_start + site_len)
                    else:

                        guides['End'].append(self.SNV_chr_pos + (search_window - target_start))
                        guides['Start'].append(self.SNV_chr_pos + (search_window - target_start) - site_len)

                    guides['Chr'].append(self.chrom)
                    guides['Strand'].append(strand)
                    guides['Pam'].append(str(guide[20:]))
                    guides['gRNA'].append(str(guide[0:20]))
                    guides['SNV Position'].append(
                        (search_window - target_start) if strand == "-" else (search_window - target_start) - 1)

                search_seq = search_seq[found + 1:]
                slider = (2 * search_window) - len(search_seq)
                found = search_seq.find('GG')
        self.guides = guides

        return self.guides

'''
 test

dh = DataHandler()
clindf = dh.get_ClinVartable(hgvs_id = 'NM_000152.5(GAA):c.271G>T')
guides = dh.get_Guide(pam = 'NGG', win_size = [4,8],score = True, guidelen = 20)
for k,v in guides.items():
    print(k,v)

Chr ['17', '17', '17', '17', '17', '17']
Start [80104842, 80104841, 80104840, 80104839, 80104838, 80104852]
End [80104865, 80104864, 80104863, 80104862, 80104861, 80104875]
Strand ['-', '-', '-', '-', '-', '+']
gRNA ['GGCGCAAACGAAGCGGCTGT', 'GCGCAAACGAAGCGGCTGTT', 'CGCAAACGAAGCGGCTGTTG', 'GCAAACGAAGCGGCTGTTGG', 'CAAACGAAGCGGCTGTTGGG', 'CTTCGTTTGCGCCCCTGACA']
Pam ['TGG', 'GGG', 'GGG', 'GGG', 'GGG', 'AGG']
Score [6, 0, 7, 13, 7, 36]
SNV Position [8, 7, 6, 5, 4, 4]
'''

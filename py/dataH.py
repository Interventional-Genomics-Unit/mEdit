import pandas as pd
import bbi
import subprocess
from datetime import datetime, timedelta, date
import regex as re


class DataHandler:
    '''
    To call tables and handle clinvar updates
    '''

    def __init__(self):
        #paths
        self.datadir = "/groups/clinical/projects/editability/"
        self.clinvar_folder = f"{self.datadir}clinvar/"
        self.clinvar_summary = f"{self.datadir}clinvar/variant_summary.txt.gz"
        self.CRISPRpath = f"{self.datadir}CRISPR_hg38/"
        self.HGVSlookup_path = f"{self.clinvar_folder}HGVSlookup.csv"
        self.lastUpdate_file = f"{self.datadir}clinvar/clinvar_lastUpdate.txt"

        #tables
        self.guide_tab = None
        self.vardf = None

        #Variables
        self.chrom = str()
        self.hgvs_id = str()
        self.ref_allele = str()
        self.alt_allele = str()



        #Determine if an update is needed
        f = open(self.lastUpdate_file,"r").readlines()
        lastdate = datetime.strptime(f[0], '%Y-%m-%d').date()
        if (date.today() - lastdate).days > 31:
            self.updateTables()


    def updateTables(self):
        #add user inquiry to whether they want ot init update

        print("updating to latest version of clinvar")
        clinvar_ftp = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        cmd = f"wget {clinvar_ftp} -O {self.clinvar_summary}"
        p = subprocess.run(cmd, shell=True,
                           capture_output=True)
        print(p)

        with open(self.lastUpdate_file, "w") as f:
            today = date.today()
            f.write(str(today))
        f.close()
        #import Update_clinvar.py to parse newly updated filed


    def get_HGVStable(self):
        '''
        Get HGVS Lookup table - a table to find the chrom a HGVS ID resides
        '''
        hgvs_tab = pd.read_csv(self.HGVSlookup_path)
        return hgvs_tab

    def getChrom(self, hgvs_id):

        self.hgvs_id = hgvs_id
        hgvs_tab = self.get_HGVStable()
        hgvs_list = hgvs_tab["HGVS"]
        try:
            chrom = hgvs_tab['Chr'].iloc[list(hgvs_list).index(hgvs_id.split(".")[0])]
            self.chrom = chrom

        except IndexError:
            raise "No HGVSs were found for this search"
        return chrom

    def get_ClinVartable(self,prefix):
        # df = pd.read_csv("/groups/clinical/projects/editability/clinvar/17_variant.txt", dtype = object)
        df = pd.read_csv(f"{self.clinvar_folder}{self.chrom}_variant.txt", dtype = object)
        df = df.loc[df["Name"].str.startswith(prefix)]
        df1 = df[df.Name == self.hgvs_id]
        if df1.shape[0] < 1:
            try:
                var_desc = re.search("(\d+|c|.)[.+-]\d+\w+>\w+", self.hgvs_id).captures()[0]
                try:
                    df = df.loc[df["Name"].str.contains(var_desc)]
                except IndexError:
                    raise "No HGVSs were found for this search"
            except:
                print("Be sure hgvs id contains the variant description ex: 271G>A")
        else:
            df = df1

        self.vardf = df
        self.ref_allele, self.alt_allele = df.ReferenceAlleleVCF.iloc[0], df.AlternateAlleleVCF.iloc[0]
        return df



    def get_Guidetable(self, chrom, pos, win_size):
        bbfile = f"{self.CRISPRpath}chr{str(chrom)}.bb"
        #bbfile = "/groups/clinical/projects/editability/CRISPR_hg38/chr17.bb"
        #pos = 80104857
        f = bbi.open(bbfile)
        pam_len = 3
        site_len = 20 + pam_len

        #search a broader range first
        df = f.fetch_intervals(chrom=f"chr{chrom}", start=pos -site_len, end=pos + site_len)
        df = df.rename(columns={'field8': 'gRNA', 'field9': 'PAM', 'field10': 'Scores'})

        #narrow search by strand loc
        ## negative strand search
        neg_max, neg_min = (pos+1) + (win_size[1] - site_len), (pos+1) + (win_size[0] - site_len)
        df1 = df[df.strand == "-"]
        df1 = df1[df1.start >= neg_min]
        df1 = df1[df1.start < neg_max]
        if df1.shape[0] > 0:
            df1['SNV_position']= [(x - pos) for x in df1.end]
            #SWAP gRNA Ref allele to Alt
            altrev = 'G' if self.alt_allele == 'C' else 'T' if self.alt_allele == 'A' else 'A' if self.alt_allele == 'T' else 'C'
            df1['gRNA'] =df1[['gRNA', 'SNV_position']].apply(lambda x: (x.gRNA[0:x.SNV_position] + altrev + x.gRNA[x.SNV_position + 1:]),
                                               axis=1)

        ## postive strand search
        pos_min, pos_max = (pos+1) - win_size[1], (pos+1) - win_size[0]
        df2 = df[df.strand == "+"]
        df2 = df2[df2.start >=pos_min]
        df2 = df2[df2.start < pos_max]
        if df2.shape[0] > 0:
            # SWAP gRNA Ref allele to Alt
            df2['SNV_position'] = [((pos - x)-1) for x in df2.start]
            df2['gRNA'] = df2[['gRNA', 'SNV_position']].apply(lambda x: (x.gRNA[0:x.SNV_position] + self.alt_allele + x.gRNA[x.SNV_position + 1:]),
                                                axis=1)

        df = pd.concat([df1, df2], ignore_index=True).drop_duplicates().drop(columns = ['name', 'score'])

        self.guide_tab = df
        return self.guide_tab



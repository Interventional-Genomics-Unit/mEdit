################
#   Keeps Clinvar data up-to-date
#   Checks User input and Terms
###############

import pandas as pd
from datetime import datetime, date
import subprocess
import gzip
import regex as re


class Validator:
    '''
    Validates inputs and datbases
    '''

    clinvar_ftp = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    clinvar_vcf = ' https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz'
    clinvar_index = ' https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi'
    refseq = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"

    def __init__(self, datadir):

        # database paths
        self.raw_tables = f"{datadir}raw_tables/"
        self.clinvar_summary = f"{self.raw_tables}clinvar/variant_summary.txt.gz"  # raw clinvar table
        self.HPApath = f"{self.raw_tables}HPA/proteinatlas.tsv"
        self.gencode_path = f"{self.raw_tables}gencode/GENEonly_cleaned_genecode_annotation.csv"  # Genecode Table

        self.processed_tables = f"{datadir}/processed_tables/"
        self.HGVSlookup_path = f"{self.processed_tables}HGVSlookup.csv" #Chrom to refID key table
        self.lastUpdate_file = f"{self.processed_tables}clinvar_lastUpdate.txt"

        #variables
        self.chroms = ['1', '2','3', '4', '5', '6', '7', '8', '9', '10', '11','12',
                         '13', '14', '15', '16', '17', '18', '19', '20', '21', '22','MT', 'Y', 'X']
        self.possible_queryTypes = ["coordinates","hgvs","phenotype"]
        self.possible_editors = ["all", "spCas9","baseeditor"]
    def process_refseq(self):
        """
        ftp : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
        """
        mane = pd.read_csv(f'{self.processed_tables}MANE.GRch38.summary_cleaned.txt.gz')
        mane_ensid = list(mane['Ensembl_TranscriptID'])
        mane_tid = list(mane['TranscriptID'])

        labels = ['bin', 'id', 'chrom', 'strand', 'txStart', 'txEnd',
                  'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
                  'score', 'name', 'cdsStartStat', 'cdsEndStat',
                  'exonFrames']
        out = gzip.open(f'{self.processed_tables}ncbiRefSeq.txt.gz', 'wt')
        for line in gzip.open(f'{self.raw_tables}Refseq/ncbiRefSeq.txt.gz', 'rt'):
            tokens = line.split('\t')
            if tokens[2].replace('chr', "") in self.chroms:
                try:
                    i = mane_tid.index(tokens[1])
                    tokens[0] = mane_ensid[i]
                except:
                    tokens[0] = '-'
                out.write('\t'.join(tokens[0:8] + tokens[9:11] + [tokens[12]] + [tokens[-1]]))
        out.close()

    def process_MANE(self):
        '''
        Mane has joins ENS and REFSEQ IDs
        '''
        mane = pd.read_csv(f'{self.raw_tables}clinvar/MANE.GRch38.v1.1.summary.txt.gz', sep="\t")

        mane['#NCBI_GeneID'] = mane['#NCBI_GeneID'].str.replace('GeneID:',"")
        mane = mane[mane['MANE_status'] == 'MANE Select']
        mane = mane.rename(columns = {'#NCBI_GeneID': 'GeneID','symbol':'GeneSymbol','RefSeq_nuc':'TranscriptID',
                             'RefSeq_prot':'ProteinID','chr_start':'Start',
                             'chr_end':'End','chr_strand':'Strand','GRCh38_chr':'ChrID','Ensembl_nuc':'Ensembl_TranscriptID'})

        mane = mane.drop(columns = 'MANE_status')
        mane = mane.loc[mane.ChrID.str.startswith('NC_')]
        mane.to_csv(f'{self.processed_tables}MANE.GRch38.summary_cleaned.txt.gz',index = False, compression='gzip')
    def process_clinvarVCF(self):
        '''
        clinvarvcf has molecular consequences
        '''
        mc = pd.read_csv('/groups/clinical/projects/editability/tables/raw_tables/clinvar/clinvar.vcf'
                         , comment='#', sep='\t', names=['Chr', 'PositionVCF', '?', 'REF', 'ALT', 'x', 'x1','attributes'])
        mc['chrHGVS_ID'] = mc['attributes'].str.extract(r'CLNHGVS=([^;]*)', expand=False)
        mc['MC'] = mc['attributes'].str.extract(r'MC=([^;]*)', expand=False)
        mc['AlleleID'] = mc['attributes'].str.extract(r'ALLELEID=([^;]*)', expand=False)
        con = []
        mc.loc[mc['MC'].isna(),'MC'] = '-'
        for x in mc['MC']:
            if '|' in x:
                x = ','.join([y for y in x.split('|') if 'SO' not in y])
            con.append(x)
        mc['MC'] = con

        mc = mc.drop(columns=['x1', 'x', 'attributes'])
        mc.to_csv(f'{self.processed_tables}clinvarvcf2txt.txt', index=False)

    def add_molecular_consequences(self):
        '''
        add molecular consequences
        '''
        for ch in self.chroms:
            df = pd.read_csv(f"{self.processed_tables}variant_tables/{ch}_variant.txt")
            mc = pd.read_csv(f'{self.processed_tables}clinvarvcf2txt.txt')
            mc['Chr'] = mc['Chr'].astype('str')
            mc = mc[mc['Chr'] == ch]
            mc = mc[['MC', 'AlleleID']]
            mc['AlleleID'] = mc['AlleleID'].astype('int64')
            joined = df.join(mc.set_index('AlleleID'), on='AlleleID')
            joined.to_csv(f'{self.processed_tables}variant_tables/{ch}_variant.txt', index=False)


    def extract_tid_from_hgvs(self,vdf):
        rprefix = r"((N(M|G|C|R)_[\d.]*)|(m))"
        rsuffix = r"(:(c|m|g|n)\.\S*)"
        hgs = list(vdf['HGVS_ID'])
        simple_ids = []
        tids = []
        for h in hgs:
            if re.search(rsuffix,h) and re.search(rprefix,h):
                tid = re.search(rprefix,h).groups()[0]
                suf = tid + re.search(rsuffix,h).groups()[0]
                tids.append(tid)
                simple_ids.append(suf)
            elif h.startswith('m'):
                tid = 'm'
                suf = h
                tids.append(tid)
                simple_ids.append(suf)

            else:
                print(h)
                tids.append('-')
                simple_ids.append('-')

        vdf.insert(3, 'HGVS_Simple', simple_ids)
        vdf.insert(4,'TranscriptID',tids)
        return vdf

    def add_MANE(self):
        mane = pd.read_csv(f'{self.processed_tables}MANE.GRch38.summary_cleaned.txt.gz')
        mane = mane[['Start', 'End', 'Strand', 'TranscriptID', 'ProteinID', 'Ensembl_TranscriptID', 'Ensembl_Gene']]

        for ch in self.chroms:
            vdf = pd.read_csv(f"{self.processed_tables}variant_tables/{ch}_variant.txt")
            mane['Ensembl_Gene'] = [x.split('.')[0] if type(x) == str else x for x in mane['Ensembl_Gene']]
            joined = vdf.join(mane.set_index('TranscriptID'), on = 'TranscriptID')

            joined.to_csv(f"{self.processed_tables}variant_tables/{ch}_variant.txt", index=False)
        #colreorder = [2,4,16,29,30,31,25,26,0,1,3,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,27,28,32,33,34]
        #joined = joined.iloc[:,colreorder]


    def extractOMIM(self,vdf):
        '''
        Find OMIM ID clinvar table and make seperate column
        '''
        all_ids = [",".join([x, y]) if len(",".join([x, y])) > 0 else "NA" for x, y in
                   zip(vdf['PhenoIDS'], vdf['OtherIDs'])]
        omim = []
        new_all_ids = []
        for x in all_ids:
            x = x.replace("MONDO:MONDO:", "MONDO:")
            found = list(set(re.findall("OMIM:([PS\.0-9]+)", x)))
            if len(found) == 0:
                omim.append("-")
            elif len(found) == 1:
                omim.append(found[0])
                x = x.replace(f"OMIM:{found[0]}", "")
            else:
                omim.append("|".join([z for z in found]))
                for z in found:
                    x = x.replace(f"OMIM:{z}", "")
            new_all_ids.append(x)
        vdf["OMIM"] = omim
        vdf["IDs"] = new_all_ids
        vdf = vdf.drop(columns=['PhenoIDS', 'OtherIDs'])
        return vdf

    def clean_clinvar(self):
        '''
        splits clinvar into files by chromosomes,
        removes unneeded columns
        Keeps only data from hg38 assembly
        '''
        ## Dropped cols
        # "LastEvaluated", "RS(dbSNP)", "Origin", 'Chromosome',
        # , "Cytogenetic", "ReviewStatus",
        # "NumberSubmitters", "Gudelines", "TestedInGTR", "SubmitterCategories"
        to_drop = ['RefAllele', 'AltAllele','Assembly', 'Start', 'Stop']
        in_file = gzip.open(self.clinvar_summary, "rt")
        contents = in_file.readlines()

        allcols = ['AlleleID', 'Type', 'HGVS_ID', 'GeneID', 'GeneSymbol', 'HGNC_ID',
                   'ClinicalSign', "ClinSigSimple", "LastEval", "RS#(dbSNP)", "nsv/esv (dbVar)",
                   "RCVaccession", "PhenoIDS","PhenoList", "Origin", "OriginSimple", "Assembly",
                   "ChrAccession", "Chr", "Start","Stop","RefAllele", "AltAllele",
                   "Cytogenetic", "ReviewStatus", "NumberSubmitters", "Guidelines", "TestedInGTR", "OtherIDs",
                   "SubmitterCategories", "VariationID", "PositionVCF", "RefAlleleVCF", "AltAlleleVCF"]

        cols = [allcols.index(c) for c in allcols if c not in to_drop]
        for ch in self.chroms:
            print(ch)
            out_fname = f"{self.processed_tables}/variant_tables/{ch}_variant.txt"
            lines = []
            for line in contents:
                line = line.split("\t")
                if line[18] == str(ch):
                    if line[16] != "GRCh37": # Remove hg19 data
                        line[-1] = line[-1].replace("\n", "")
                        if 'single nucleotide variant' in line:
                            line = [line[i] for i in cols]
                            lines.append(line)

            vdf = pd.DataFrame(lines, columns = [allcols[i] for i in cols])
            vdf['HGNC_ID'] = vdf['HGNC_ID'].apply(lambda x: x.replace("HGNC:",""))
            vdf['PositionVCF'] = vdf['PositionVCF'].astype('int')
            vdf = vdf[vdf['PositionVCF']>1]

            #keep only pathogenic
            vdf = vdf[vdf['ClinSigSimple'] == '1']##variants remaining = 12426
            vdf = vdf.sort_values('GeneSymbol')

            #Extract GeneID and GeneSymbol from HGVSID, In places that are missing
            blankgenes = vdf.loc[vdf['GeneID'] == '-1']
            hgvs = blankgenes['HGVS_ID']
            genenames = hgvs.str.extract(r"\((\w*\d*)\):").sort_values(0)
            temp = vdf.loc[~vdf['HGVS_ID'].isin(blankgenes)].drop_duplicates(subset='GeneSymbol')
            genenames = genenames.join(temp[['GeneID','GeneSymbol']].set_index('GeneSymbol'),on = 0)
            vdf.loc[vdf['HGVS_ID'].isin(hgvs),'GeneSymbol'] = genenames[0]
            vdf.loc[vdf['HGVS_ID'].isin(hgvs),'GeneID'] = genenames['GeneID']
            vdf.loc[vdf.GeneID.isna(),'GeneID'] = -1
            vdf['GeneID'] = vdf['GeneID'].astype('int64')
            vdf = vdf.reset_index(drop=True)

            # Find and extract OMIM ID
            vdf = self.extract_tid_from_hgvs(vdf)
            vdf = self.extractOMIM(vdf)
            vdf.to_csv(out_fname,index=False)

        self.add_MANE()
        self.add_molecular_consequences()


    def intervalMatch(self,df1,df2):
        '''
        matches a snp position to gene coord by creating dummy indexes
        '''
        df1 = df1.reset_index() #clinvar
        df2 = df2.reset_index()  # hpa
        temp1 = [f'unmatched{str(x)}' for x in df1["PositionVCF"].index] #dummy clin var index
        temp2 = [f'matched{x}' for x in df2["Position"].index] #dummy hpa index
        starts = [int(x.split("-")[0]) for x in df2.Position]
        ends = [x.split("-")[1] for x in df2.Position]
        for i in df1.index:
            posvcf = df1["PositionVCF"].iloc[i]
            new_end = [y for x, y in zip(starts,ends) if int(posvcf) > int(x)]
            new_p2 = [y for y in new_end if int(posvcf) < int(y)]
            if len(new_p2)>0:
                # set clin var dummy interval to match hpa
                temp1[i] = f"matched{ends.index(new_p2[0])}"
        return temp1,temp2

    def make_gene_tables(self):
        '''
        attached HPA gene expression info to clinvar info by ensembl id
        '''

        hpa = pd.read_csv(self.HPApath, delimiter='\t')

        cols = ['Ensembl', 'Chromosome', 'Gene description', 'Protein class', 'Biological process',
                'Molecular function', 'Uniprot',
                'Disease involvement', 'RNA tissue specificity', 'RNA tissue specific nTPM',
                'RNA tissue distribution', 'RNA tissue cell type enrichment',
                'RNA single cell type specific nTPM', 'RNA tissue specific nTPM']
        hpa = hpa[cols]

        mane = pd.read_csv(f'{self.processed_tables}MANE.GRch38.summary_cleaned.txt.gz')
        mane = mane[['TranscriptID','Start', 'End', 'Strand',  'ProteinID', 'Ensembl_TranscriptID', 'Ensembl_Gene']]
        mane['Ensembl_Gene'] = [x.split('.')[0] if type(x) == str else x for x in mane['Ensembl_Gene']]
        joined_df = hpa.join(mane.set_index('Ensembl_Gene'), on='Ensembl_Gene')
        joined_df.to_csv(f'{self.processed_tables}gene_tables.csv.gz', index=False,compression='gzip')

    def make_HGVStable(self):
        '''
        Create CSV file of unduplicated HGVS prefix(coding ref name) and Chromosome
        '''
        names, chrs = [], []
        for ch in self.chroms:
            clin = pd.read_csv(f"{self.processed_tables}variant_tables/{ch}_variant.txt")
            names += list(set(clin['TranscriptID']))
            chrs += [ch for i in range(len(set(clin['TranscriptID'])))]

        df = pd.DataFrame({"TranscriptID": names, "Chr": chrs}).drop_duplicates()
        out_file = f"{self.processed_tables}HGVSlookup.csv"
        df.to_csv(out_file, index=None)

    def updateTables(self):
        # TODO: add user inquiry to whether they want to init update
        print("updating to latest version of clinvar")

        cmd = f"wget {self.clinvar_ftp} -O {self.clinvar_summary}"
        p = subprocess.run(cmd, shell=True,
                           capture_output=True)
        print(p)

        cmd = f"wget {self.refseq} -O {self.raw_tables}Refseq/ncbiRefSeq.txt.gz"
        p = subprocess.run(cmd, shell=True,
                           capture_output=True)
        print(p)

        cmd = f"wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.1.summary.txt.gz -O " \
              f"{self.raw_tables}clinvar/MANE.GRch38.v1.1.summary.txt.gz"
        p = subprocess.run(cmd, shell=True,
                           capture_output=True)
        print(p)

        cmd = f"wget {self.clinvar_vcf} -O {self.raw_tables}clinvar/clinvar.vcf.gz"
        p = subprocess.run(cmd, shell=True,
                           capture_output=True)
        print(p)
        cmd = f"wget {self.clinvar_index} -O {self.raw_tables}clinvar/clinvar.vcf.gz.tbi"
        p = subprocess.run(cmd, shell=True,
                           capture_output=True)

        cmd = f"bcftools view -f type!=snp {self.clinvar_vcf} -o {self.raw_tables}clinvar/clinvar.vcf -O v"
        p = subprocess.run(cmd, shell=True,
                           capture_output=True)

        print("Cleaning and Splitting Clinvar.....")
        self.clean_clinvar()
        print("Adding Ensembl Identifiers....")
        self.add_MANE()
        print("Adding Molecular Consequences....")
        self.add_molecular_consequences()
        print("Appending HPA data.....")
        self.appendHPA()
        print("Writing new HGVS Lookup table.....")
        self.make_HGVStable()

        with open(self.lastUpdate_file, "w") as f:
            today = date.today()
            f.write(str(today))
        f.close()

#     def check_updates(self):
#         """
#         Determines if an update is needed based on checking the date of txt file
#         '''
#
#         f = open(self.lastUpdate_file, "r").readlines()
#         lastdate = datetime.strptime(f[0], '%Y-%m-%d').date()
#         if (date.today() - lastdate).days > 31:
#             self.updateTables()
#         else:
#             print('You are using the latest clinvar data')
#             print(f"Last updated {lastdate}")
#
#     def validate_input_file(self, input_file):
#         cols = ["Query", "Type","Editor"]
#         try:
#             self.input_df = pd.read_csv(input_file,usecols= cols)
#         except:
#             raise FileNotFoundError(f"There was was a problem reading {input_file}. \n"
#                                     f"Be sure this is a csv file with the column headers {cols}")
#
#         for term in self.input_df['Type'].unique():
#             if term not in self.possible_queryTypes:
#                 raise ValueError(f"Query type must be either {self.possible_queryTypes}")
#
#         for term in self.input_df['Editor'].unique():
#             if term not in self.possible_editors:
#                 raise ValueError(f"Editor type must be either {self.possible_editors}")
#         return self.input_df
#
#
# '''
# datadir = "/groups/clinical/projects/editability/tables/"
# val = Validator(datadir)
# val.clean_clinvar()
# val.add_MC(
#
# val.make_guidetabs()
#
# for ch in val.chroms:
#     vdf = pd.read_csv(f"{val.processed_tables}{ch}_variant.txt")
# '''
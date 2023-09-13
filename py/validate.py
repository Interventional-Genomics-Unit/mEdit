################
#   Keeps Clinvar data up-to-date
#   Checks User input and Terms
###############
import pandas as pd
from datetime import datetime, date
import subprocess
import gzip
import regex as re


#datadir = "/groups/clinical/projects/editability/tables/"

class Validator:
    '''
    Validates inputs and datbases
    '''

    clinvar_ftp = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"


    def __init__(self, datadir):

        # database paths
        self.raw_tables = f"{datadir}/raw_tables/"
        self.clinvar_summary = f"{self.raw_tables}variant_summary.txt.gz"  # raw clinvar table
        self.gencode_path = f"{self.raw_tables}gencode/GENEonly_cleaned_genecode_annotation.csv" #Genecode Table
        self.HPApath = f"{self.raw_tables}HPA/proteinatlas.tsv"

        self.processed_tables = f"{datadir}/processed_tables/"
        self.HGVSlookup_path = f"{self.processed_tables}HGVSlookup.csv" #Chrom to refID key table
        self.lastUpdate_file = f"{self.processed_tables}clinvar_lastUpdate.txt"


        #user input
        self.input_df = pd.DataFrame()

        #variables
        self.chroms = [ '18', '3', 'MT', 'Y', 'X']
        self.possible_queryTypes = ["coordinates","hgvs","phenotype"]
        self.possible_editors = ["all", "spCas9","baseeditor"]

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

    def add_Gencode(self,vdf,ch):

        gdf = pd.read_csv(self.gencode_path)
        gdf["GeneSymbol"] = gdf['Gene_Name']
        gdf = gdf[gdf['Chr'] == str(ch)]
        gdf['End'] = gdf['End'].astype('int64')
        gdf['Start'] = gdf['Start'].astype('int64')
        gdf = gdf.drop_duplicates(subset = 'GeneSymbol')
        vdf['GeneSymbol'] = vdf['GeneSymbol'].str.replace("[;]?LOC\d*[;]?", "", regex = True) #remove non-specific secondary names
        not_matching_gencode = list(set(vdf.GeneSymbol).difference(set(gdf.GeneSymbol)))
        matched = vdf.loc[~vdf["GeneSymbol"].isin(not_matching_gencode)]
        joined_df = matched.join(gdf.set_index('GeneSymbol'), on='GeneSymbol').reset_index(drop = True)

        remainder = vdf.loc[vdf["GeneSymbol"].isin(not_matching_gencode)].reset_index(drop = True)

        new_rows = []
        for i in remainder.index:
            posvcf = remainder["PositionVCF"].iloc[i]
            row = gdf[gdf['Start'].lt(posvcf) & gdf['End'].gt(posvcf)]

            if len(row['Gene_Name']) == 0:
                new_rows.append(list(remainder.iloc[i]) +[ch,int(),int(),"","","",""])

            if len(row['Gene_Name']) > 0:
                for x in range(len(row.index)):
                    if 'ENS' not in row.iloc[x,4]:
                        new_rows.append(list(remainder.iloc[i]) + list(row.iloc[x,0:7]))

        joined_df = pd.concat([joined_df, pd.DataFrame(new_rows, columns = joined_df.columns)])
        joined_df = joined_df.iloc[:,[2,4,19,20,21,22,16,0,1,3,5,6,7,8,9,10,11,12,13,14,15,17,18,23,24,25]]
        return joined_df

    def clean_clinvar(self):
        '''
        splits clinvar into files by chromosomes,
        removes unneeded columns
        Keeps only data from hg38 assembly
        '''
        ## Dropped cols
        # "LastEvaluated", "RS(dbSNP)", "Origin", 'Assembly','Chromosome','Start', 'Stop',
        # 'ReferenceAllele', 'AlternateAllele', "Cytogenetic", "ReviewStatus",
        # "NumberSubmitters", "Gudelines", "TestedInGTR", "SubmitterCategories"
        to_drop = [8, 9, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29]
        in_file = gzip.open(self.clinvar_summary, "rt")
        contents = in_file.readlines()

        allcols = ['AlleleID', 'Type', 'HGVS_ID', 'GeneID', 'GeneSymbol', 'HGNC_ID',
                   'ClinicalSign', "ClinSigSimple", "LastEval", "RS#(dbSNP)", "nsv/esv (dbVar)",
                   "RCVaccession", "PhenoIDS","PhenoList", "Origin", "OriginSimple", "Assembly",
                   "ChrAccession", "Chr", "Start","Stop","RefAllele", "AltAllele",
                   "Cytogenetic", "ReviewStatus", "NumberSubmitters", "Guidelines", "TestedInGTR", "OtherIDs",
                   "SubmitterCategories", "VariationID", "PositionVCF", "RefAlleleVCF", "AltAlleleVCF"]

        cols = [allcols[i] for i in range(34) if i not in to_drop]
        for ch in self.chroms:
            out_fname = f"{self.processed_tables}{ch}_variant.txt"
            lines = []
            for line in contents:
                line = line.split("\t")
                if line[18] == str(ch):
                    if line[16] != "GRCh37": # Remove hg19 data
                        line[-1] = line[-1].replace("\n", "")
                        if len(line[-1]) < 2 and line[-1] != 'na': # remove Alt allele that are less than 2bp
                            if len(line[-2]) < 2 and line[-2] != 'na':  # remove Alt allele that are less than 2bp
                                line = [line[i] for i in range(34) if i not in to_drop]
                                lines.append(line)

            vdf = pd.DataFrame(lines, columns = cols)
            vdf['HGNC_ID'] = vdf['HGNC_ID'].apply(lambda x: x.replace("HGNC:",""))
            vdf['PositionVCF'] = vdf['PositionVCF'].astype('int').sort_values()
            vdf = vdf[vdf['PositionVCF']>1]

            #Find and extract OMIM ID
            vdf = self.extractOMIM(vdf)

            #annotate gene info with genecode
            vdf = self.add_Gencode(vdf,ch)

            vdf.to_csv(out_fname,index=None)


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


    def appendHPA(self):
        '''
        attached HPA gene expression info to clinvar info by using gene names first and then by chromosome locations
        '''

        hpa_og = pd.read_csv(self.HPApath, delimiter='\t')

        cols = ['Gene', 'Gene synonym', 'Ensembl', 'Gene description',
       'Chromosome', 'Position', 'Biological process', 'Molecular function',
       'Disease involvement', 'RNA tissue specific nTPM',
       'RNA tissue distribution', 'RNA tissue cell type enrichment',
       'RNA single cell type specific nTPM', 'RNA tissue specific nTPM']
        hpa_og = hpa_og[cols]

        for ch in self.chroms:
            print(ch)
            #import cleaned clinvar
            vdf = pd.read_csv(f"{self.processed_tables}{ch}_variant.txt")
            hpa = hpa_og[hpa_og["Chromosome"] == str(ch)]
            vdf.Gene_ID = [str(x).split(".")[0] for x in vdf.Gene_ID]
            vdf = vdf.rename(columns = {'Gene_ID': 'Ensembl',"strand": "Coding Strand"})

            #Try to match on ensembl ID first
            ensembl_diff= list(set(vdf.Ensembl).difference(set(hpa.Ensembl)))
            not_e_match = vdf.loc[vdf['Ensembl'].isin(ensembl_diff)]
            e_match = vdf.loc[~vdf['Ensembl'].isin(ensembl_diff)]
            joined_df = e_match.join(hpa.set_index('Ensembl'),on='Ensembl')

            #for the remainder of unmatched find genes by coords
            #This is too slow to run for the entire dataset
            temp1,temp2 = self.intervalMatch(df1=not_e_match, df2=hpa) #creates dummy matching indexes
            not_e_match['matched'] = temp1
            hpa['matched'] = temp2
            joined_df2 = not_e_match.join(hpa.set_index('matched').drop(columns = 'Ensembl'), on='matched')
            df = pd.concat([joined_df, joined_df2.drop(columns=['matched'])], ignore_index=True).reset_index(drop = True)
            df['GeneSymbol'] = df['GeneSymbol'].str.replace(";|:|\\|", ",", regex=True)
            df['Gene_Name'] = df['Gene_Name'].str.replace(";|:|\\|", ",", regex=True)
            df['Gene synonym'] = df['Gene synonym'].str.replace(";|:|\\|", ",", regex=True)
            df['Gene'] = df['Gene'].str.replace(";|:|\\|", ",", regex=True)
            gene_syns = []
            for i in df.index:
                names = f"{df['GeneSymbol'].iloc[i]},{str(df['Gene_Name'].iloc[i])},{str(df['Gene synonym'].iloc[i])}" \
                        f",{str(df['Gene'].iloc[i])}"
                names = set(names.replace("nan,","").strip().split(","))
                gene_syns.append(";".join(i for i in names))
            df['Gene synonym'] = gene_syns
            df = df.drop(columns =['Gene_Name','Gene'])

            df.to_csv(f'{self.processed_tables}{ch}_variant.txt', index=False)
        return df



    def make_HGVStable(self):
        '''
        Create CSV file of unduplicated HGVS prefix(coding ref name) and Chromsome
        '''
        names, chrs = [], []
        for ch in self.chroms:
            clin = pd.read_csv(f"{self.processed_tables}{ch}_variant.txt")
            names += list(clin['HGVS_ID'])
            chrs += [ch for i in range(len(clin['HGVS_ID']))]

        df = pd.DataFrame({"HGVS": names, "Chr": chrs})
        new_HGVS = [x.split(".")[0] for x in df['HGVS']]
        df = pd.DataFrame({"HGVS": new_HGVS, "Chr": chrs}).drop_duplicates()
        df = df.loc[df["HGVS"].str.startswith("NM")]
        out_file = f"{self.processed_tables}HGVSlookup.csv"
        df.to_csv(out_file, index=None)



    def updateTables(self):
        # TODO: add user inquiry to whether they want to init update
        print("updating to latest version of clinvar")

        cmd = f"wget {self.clinvar_ftp} -O {self.clinvar_summary}"
        p = subprocess.run(cmd, shell=True,
                           capture_output=True)
        print(p)

        print("Cleaning and Splitting Clinvar.....")
        self.clean_clinvar()
        print("Appending HPA data.....")
        self.appendHPA()
        print("Writing new HGVS Lookup table.....")
        self.make_HGVStable()


        with open(self.lastUpdate_file, "w") as f:
            today = date.today()
            f.write(str(today))
        f.close()


    def check_updates(self):
        '''
        Determines if an update is needed based on checking the date of txt file
        '''

        f = open(self.lastUpdate_file, "r").readlines()
        lastdate = datetime.strptime(f[0], '%Y-%m-%d').date()
        if (date.today() - lastdate).days > 31:
            self.updateTables()
        else:
            print('You are using the latest clinvar data')
            print(f"Last updated {lastdate}")

    def validate_input_file(self, input_file):
        cols = ["Query", "Type","Editor"]
        try:
            self.input_df = pd.read_csv(input_file,usecols= cols)
        except:
            raise FileNotFoundError(f"There was was a problem reading {input_file}. \n"
                                    f"Be sure this is a csv file with the column headers {cols}")

        for term in self.input_df['Type'].unique():
            if term not in self.possible_queryTypes:
                raise ValueError(f"Query type must be either {self.possible_queryTypes}")

        for term in self.input_df['Editor'].unique():
            if term not in self.possible_editors:
                raise ValueError(f"Editor type must be either {self.possible_editors}")


        return self.input_df

    def validate_query_terms(self):
        #TODO Validated HGVS terms or validate
        pass



import pandas as pd
import gzip
import regex as re
import os
import subprocess


#----------------------------------------#
# uploads newest clinvar files
# sorts and cleans clinvar
# appends HPA file info
#----------------------------------------#


def split_clinVar(clinvar_file):
    '''
    splits clinvar into files by chromosomes,
    removes unneeded columns
    Keeps only data from hg38 assembly
    Elimnates
    '''
    chroms = ['7','15', '11', '14', '6', '2', '20', '10', '19', '16', '22', '12',
              '1', '8', '9', '13', '21', '5', '4', '17', '18', '3', 'MT','Y', 'X']
    to_drop = [8,9,14,16,18,19,20,21,22,23,24,25,26,27,29]

    ## Dropped cols = "LastEvaluated", "RS(dbSNP)", "Origin", 'Assembly','Chromosome','Start', 'Stop', 'ReferenceAllele', 'AlternateAllele',
    # "Cytogenetic", "ReviewStatus", "NumberSubmitters",
    ##"Gudelines", "TestedInGTR", "SubmitterCategories"]
    in_file = gzip.open(clinvar_file, "rt")
    contents = in_file.readlines()
    allcols = ['AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol', 'HGNC_ID',
               'ClinicalSignificance', "ClinSigSimple", "LastEvaluated", "RS# (dbSNP)", "nsv/esv (dbVar)",
               "RCVaccession", "PhenotypeIDS","PhenotypeList", "Origin", "OriginSimple", "Assembly",
               "ChromosomeAccession", "Chromosome", "Start","Stop","ReferenceAllele", "AlternateAllele",
               "Cytogenetic", "ReviewStatus", "NumberSubmitters", "Gudelines", "TestedInGTR", "OtherIDs",
               "SubmitterCategories", "VariationID", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
    cols = [allcols[i] for i in range(34) if i not in to_drop]
    path = clinvar_file[:clinvar_file.rfind["/"]+1]
    for ch in chroms:
        out_fname = f"{path}{ch}_variant.txt"
        lines = []

        for line in contents:
            line = line.split("\t")
            if line[18] == str(ch):
                if line[16] != "GRCh37": # Remove hg19 data
                    line[-1] = line[-1].replace("\n", "")
                    if len(line[-1]) < 2 and line[-1] != 'na': # remove Alt allele that are less than 2bp
                        line = [line[i] for i in range(34) if i not in to_drop]
                        lines.append(line)

        vdf = pd.DataFrame(lines,columns = cols)
        vdf = vdf[vdf.Type != 'copy number gain']
        vdf = vdf[vdf.Type != 'copy number loss']
        vdf['HGNC_ID'] = vdf['HGNC_ID'].apply(lambda x: x.replace("HGNC:",""))
        #Find OMIM ID
        all_ids = [",".join([x,y]) if len(",".join([x,y]) ) > 0 else "NA" for x,y in zip(vdf['PhenotypeIDS'],vdf['OtherIDs'])]
        omim = []
        new_all_ids = []
        for x in all_ids:
            x = x.replace("MONDO:MONDO:","MONDO:")
            found = list(set(re.findall("OMIM:([PS\.0-9]+)", x)))
            if len(found) == 0:
                omim.append("-")
            elif len(found) == 1:
                omim.append(found[0])
                x = x.replace(f"OMIM:{found[0]}","")
            else:
                omim.append("|".join([z for z in found]))
                for z in found:
                    x = x.replace(f"OMIM:{z}","")
            new_all_ids.append(x)
        vdf["OMIM"] = omim
        vdf["IDs"] = new_all_ids
        vdf = vdf.drop(columns = ['PhenotypeIDS','OtherIDs'])
        vdf['PositionVCF'] = vdf['PositionVCF'].astype('int')
        vdf = vdf.sort_values(['PositionVCF'])
        vdf = vdf[vdf['PositionVCF']>1]
        vdf.to_csv(out_fname,index=None)
    return vdf


def intervalMatch(df1,df2):
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



def appendHPA(hpafile, clinvar_folder = "/groups/clinical/projects/editability/clinvar/"):
    '''
    attached HPA gene expression info to clinvar info by using gene names first and then by chromosone locations
    '''

    hpa_og = pd.read_csv(hpafile,delimiter='\t')
    cols = ['Gene', 'Gene synonym', 'Ensembl', 'Gene description', 'Uniprot',
            'Chromosome', 'Position', 'Biological process',
            'Molecular function', 'Disease involvement',
            'RNA tissue distribution',
            'RNA tissue cell type enrichment','RNA single cell type specific nTPM',
            'RNA tissue specific nTPM', 'RNA blood cell specific nTPM']
    hpa_og = hpa_og[cols]
    chroms = ['7', '15', '11', '14', '6', '2', '20', '10', '19', '16', '22', '12',
              '1', '8', '9', '13', '21', '5', '4', '17', '18', '3', 'MT', 'Y', 'X']


    for ch in chroms:
        #import cleaned clinvar
        if clinvar_folder.endswith("/") == False:
            clinvar_folder = clinvar_folder + "/"

        clin = pd.read_csv(f"{clinvar_folder}{ch}_variant.txt")
        clin = clin[clin.GeneID != '-']

        hpa = hpa_og[hpa_og["Chromosome"] == str(ch)]

        #extract clinvar gene names that do not have a match in HPA
        clindiff = list(set(clin.GeneSymbol).difference(set(hpa.Gene)))

        # Create dataframe for only those expected to have a match
        unmatching_i = clin.loc[clin["GeneSymbol"].isin(clindiff)].index
        df_matched = clin.drop(unmatching_i)

        # Create dataframe for that genes that will not match
        matching_i = clin.loc[~clin["GeneSymbol"].isin(clindiff)].index
        df_unmatched = clin.drop(matching_i)

        # join all matching by gene symbol/name
        joined_df = df_matched.join(hpa.set_index('Gene'), on = 'GeneSymbol')

        #for the remainder of unmatched find genes by coords
        #This is too slow to run for the entire dataset
        temp1,temp2 = intervalMatch(df1=df_unmatched, df2=hpa) #creates dummy matching indexes
        df_unmatched['matched'] = temp1
        hpa['matched'] = temp2
        joined_df2 = df_unmatched.join(hpa.set_index('matched'), on='matched')
        joined_df2 = joined_df2.drop(columns=['matched','Gene'])
        df = pd.concat([joined_df,joined_df2], ignore_index=True)
        df.to_csv(f"/groups/clinical/projects/editability/clinvar/{ch}_variant.txt", index=None)



def updateClinVarHPA(clinvar_ftp,out_dir = "/groups/clinical/projects/editability/clinvar",HPA_file="/groups/clinical/projects/editability/HPA_RNA_TissueExp/proteinatlas.tsv"):
    # upload newest clinvar
    if out_dir.endswith("/") == False:
        out_dir = out_dir + "/"
    file = out_dir + "variant_summary.txt.gz"
    cmd = f"wget {clinvar_ftp} -O {file}"
    p = subprocess.run(cmd, shell=True,
                       capture_output=True)
    print(p.stderr)

    #split clinvar file and cleanup
    split_clinVar(file)

    #append HPA data
    appendHPA(HPA_file,out_dir)




updateClinVarHPA(clinvar_ftp = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
                 out_dir = "/groups/clinical/projects/editability/clinvar/",
                 HPA_file="/groups/clinical/projects/editability/HPA_RNA_TissueExp/proteinatlas.tsv")








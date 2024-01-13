# Native modules
from datetime import datetime, date
import subprocess
import gzip
import regex as re
# Installed Modules
import pandas as pd
import yaml

# print(f"Loading configuration file {config_path}")
# 	with open(config_path, "r") as f:
# 		config = yaml.load(f, Loader=yaml.FullLoader)

'''
clinvar_ftp = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
clinvar_vcf = ' https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz'
clinvar_index = ' https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi'
refseq = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"
'''
datadir = "/groups/clinical/projects/editability/tables/"

# database paths
raw_tables = f"{datadir}raw_tables/"
clinvar_summary = f"{raw_tables}clinvar/variant_summary.txt.gz"  # raw clinvar table
mane_path = f'{raw_tables}clinvar/MANE.GRch38.v1.1.summary.txt.gz'
HPApath = f"{raw_tables}HPA/proteinatlas.tsv"


processed_tables = f"{datadir}processed_tables/"
simple_tables = f"{processed_tables}guide_acquisition_tables/"
HGVSlookup_path = f"{processed_tables}HGVSlookup.csv" #Chrom to refID key table
lastUpdate_file = f"{processed_tables}clinvar_lastUpdate.txt"
mane_cleaned = f'{processed_tables}MANE.GRch38.summary_cleaned.txt.gz'



#variables
chroms = ['1', '2','3', '4', '5', '6', '7', '8', '9', '10', '11','12',
                 '13', '14', '15', '16', '17', '18', '19', '20', '21', '22','MT', 'Y', 'X']


def process_refseq(raw_tables,processed_tables):
    '''
    ftp : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz

    Makes
    A) a bed file from refseqs = 282,614 lines
    B) a bed file of the most recent genes = 66688 lines
    '''
    #in
    mane = pd.read_csv(f'{raw_tables}clinvar/MANE.GRch38.v1.1.summary.txt.gz', sep="\t")

    #out
    ref_out = gzip.open(f'{processed_tables}ncbiRefSeq.bed.gz', 'wt')

    labels = ['bin', 'id', 'chrom', 'strand', 'txStart', 'txEnd',
              'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
              'score', 'name', 'cdsStartStat', 'cdsEndStat','exonFrames']

    mane_ensid = list(mane['Ensembl_nuc'])
    mane_tid = list(mane['RefSeq_nuc'])

    latest_refs = {}
    cnt = 0
    for line in gzip.open(f'{raw_tables}Refseq/ncbiRefSeq.txt.gz', 'rt'):
        tokens = line.split('\t')
        tid, chrom, strand,tstart,tend = tokens[1:6]
        cds_start, cds_end = tokens[6], tokens[7]
        exons_start, exon_end = tokens[9], tokens[10]
        gname, frames = tokens[12], tokens[-1].split('\n')[0]

        if tokens[2].replace('chr', "") in chroms:
            tid = tokens[1]
            cnt+=1
            if tid.startswith('X') == False:
                try:
                    i = mane_tid.index(tid)
                    eid = mane_ensid[i]
                except:

                    eid = '-'
                line_out = [chrom,tstart,tend,'|'.join([strand,tid,eid,gname,cds_start,cds_end, exons_start,exon_end,frames])]
                ref_out.write('\t'.join(line_out) + '\n')
                if cnt < 20:
                    print('\t'.join(line_out))

    ref_out.close()


def process_MANE(raw_tables,processed_tables):
    '''
    Mane has joins ENS and REFSEQ IDs
    '''
    mane = pd.read_csv(f'{raw_tables}clinvar/MANE.GRch38.v1.1.summary.txt.gz', sep="\t")

    mane['#NCBI_GeneID'] = mane['#NCBI_GeneID'].str.replace('GeneID:',"")
    mane = mane[mane['MANE_status'] == 'MANE Select']
    mane = mane.rename(columns = {'#NCBI_GeneID': 'GeneID','symbol':'GeneSymbol','RefSeq_nuc':'TranscriptID',
                         'RefSeq_prot':'ProteinID','chr_start':'Start',
                         'chr_end':'End','chr_strand':'Strand','GRCh38_chr':'ChrID','Ensembl_nuc':'Ensembl_TranscriptID'})

    mane = mane.drop(columns = 'MANE_status')
    mane = mane.loc[mane.ChrID.str.startswith('NC_')]
    mane.to_csv(f'{processed_tables}MANE.GRch38.summary_cleaned.txt.gz',index = False, compression='gzip')


def process_clinvarVCF(processed_tables):
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
    mc.to_csv(f'{processed_tables}clinvarvcf2txt.txt', index=False)


def add_molecular_consequences(chroms,processed_tables):
    '''
    add molecular consequences
    '''
    for ch in chroms[10:]:
        df = pd.read_csv(f"{processed_tables}variant_tables/{ch}_variant.txt")
        mc = pd.read_csv(f'{processed_tables}clinvarvcf2txt.txt')
        mc['Chr'] = mc['Chr'].astype('str')
        mc = mc[mc['Chr'] == ch]
        mc = mc[['MC', 'AlleleID']]
        mc['AlleleID'] = mc['AlleleID'].astype('int64')
        joined = df.join(mc.set_index('AlleleID'), on='AlleleID')
        joined.to_csv(f'{processed_tables}variant_tables/{ch}_variant.txt', index=False)


def extract_tid_from_hgvs(vdf):
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
    vdf = vdf.loc[vdf['HGVS_Simple'] !='-']
    vdf = vdf.reset_index(drop=True)
    return vdf


def add_MANE(processed_tables,chroms):
    mane = pd.read_csv(f'{processed_tables}MANE.GRch38.summary_cleaned.txt.gz')
    mane = mane[['Start', 'End', 'Strand', 'TranscriptID', 'ProteinID', 'Ensembl_TranscriptID', 'Ensembl_Gene']]

    for ch in chroms:
        vdf = pd.read_csv(f"{processed_tables}variant_tables/{ch}_variant.txt")
        #mane['Ensembl_Gene'] = [x.split('.')[0] if type(x) == str else x for x in mane['Ensembl_Gene']]
        joined = vdf.join(mane.set_index('TranscriptID'), on = 'TranscriptID')

        joined.to_csv(f"{processed_tables}variant_tables/{ch}_variant.txt", index=False)
    #colreorder = [2,4,16,29,30,31,25,26,0,1,3,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,27,28,32,33,34]
    #joined = joined.iloc[:,colreorder]

def extractOMIM(vdf):
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


def clean_clinvar(clinvar_summary,chroms):
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
    in_file = gzip.open(clinvar_summary, "rt")
    contents = in_file.readlines()

    allcols = ['AlleleID', 'Type', 'HGVS_ID', 'GeneID', 'GeneSymbol', 'HGNC_ID',
               'ClinicalSign', "ClinSigSimple", "LastEval", "RS#(dbSNP)", "nsv/esv (dbVar)",
               "RCVaccession", "PhenoIDS","PhenoList", "Origin", "OriginSimple", "Assembly",
               "ChrAccession", "Chr", "Start","Stop","RefAllele", "AltAllele",
               "Cytogenetic", "ReviewStatus", "NumberSubmitters", "Guidelines", "TestedInGTR", "OtherIDs",
               "SubmitterCategories", "VariationID", "PositionVCF", "RefAlleleVCF", "AltAlleleVCF"]

    cols = [allcols.index(c) for c in allcols if c not in to_drop]
    for ch in chroms:
        print(ch)
        out_fname = f"{processed_tables}/variant_tables/{ch}_variant.txt"
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
        vdf = extract_tid_from_hgvs(vdf)
        vdf = extractOMIM(vdf)
        vdf.to_csv(out_fname,index=False)


def make_var_table(clinvar_summary,chroms):
    clean_clinvar(clinvar_summary, chroms)
    add_MANE(processed_tables, chroms)
    add_molecular_consequences(chroms,processed_tables)


def make_gene_tables(processed_tables,HPApath):
    '''
    attached HPA gene expression info to clinvar info by ensembl id
    '''

    hpa = pd.read_csv(HPApath, delimiter='\t')

    cols = ['Ensembl', 'Chromosome', 'Gene description', 'Protein class', 'Biological process',
            'Molecular function', 'Uniprot',
            'Disease involvement', 'RNA tissue specificity', 'RNA tissue specific nTPM',
            'RNA tissue distribution', 'RNA tissue cell type enrichment',
            'RNA single cell type specific nTPM', 'RNA tissue specific nTPM']
    hpa = hpa[cols]
    hpa = hpa.rename(columns = {'Ensembl':'Ensembl_Gene'})
    mane = pd.read_csv(f'{processed_tables}MANE.GRch38.summary_cleaned.txt.gz')
    mane = mane[['TranscriptID','Start', 'End', 'Strand',  'ProteinID', 'Ensembl_TranscriptID', 'Ensembl_Gene']]
    mane['Ensembl_Gene'] = [str(x).split('.')[0] for x in mane['Ensembl_Gene']]
    joined_df = hpa.join(mane.set_index('Ensembl_Gene'), on='Ensembl_Gene')
    joined_df.to_csv(f'{processed_tables}/gene_tables/gene_tables.csv.gz', index=False,compression='gzip')

def make_HGVStable(processed_tables,chroms):
    '''
    Create CSV file of unduplicated HGVS prefix(coding ref name) and Chromosome
    '''
    names, chrs = [], []
    for ch in chroms:
        clin = pd.read_csv(f"{processed_tables}variant_tables/{ch}_variant.txt")
        names += list(set(clin['TranscriptID']))
        chrs += [ch for i in range(len(set(clin['TranscriptID'])))]

    df = pd.DataFrame({"TranscriptID": names, "Chr": chrs}).drop_duplicates()
    out_file = f"{processed_tables}HGVSlookup.csv"
    df.to_csv(out_file, index=None)


def updateTables(clinvar_ftp,clinvar_summary,refseq):
    # TODO: add user inquiry to whether they want to init update
    print("updating to latest version of clinvar")

    cmd = f"wget {clinvar_ftp} -O {clinvar_summary}"
    p = subprocess.run(cmd, shell=True,
                       capture_output=True)
    print(p)

    cmd = f"wget {refseq} -O {raw_tables}Refseq/ncbiRefSeq.txt.gz"
    p = subprocess.run(cmd, shell=True,
                       capture_output=True)
    print(p)

    cmd = f"wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.1.summary.txt.gz -O " \
          f"{raw_tables}clinvar/MANE.GRch38.v1.1.summary.txt.gz"
    p = subprocess.run(cmd, shell=True,
                       capture_output=True)
    print(p)

    # cmd = f"wget {clinvar_vcf} -O {raw_tables}clinvar/clinvar.vcf.gz"
    # p = subprocess.run(cmd, shell=True,
    #                    capture_output=True)
    # print(p)
    # cmd = f"wget {clinvar_index} -O {raw_tables}clinvar/clinvar.vcf.gz.tbi"
    # p = subprocess.run(cmd, shell=True,
    #                    capture_output=True)
    #
    # cmd = f"bcftools view -f type!=snp {clinvar_vcf} -o {raw_tables}clinvar/clinvar.vcf -O v"
    # p = subprocess.run(cmd, shell=True,
    #                    capture_output=True)

    print("Cleaning and Splitting Clinvar.....")
    clean_clinvar(clinvar_summary, chroms)
    print("Adding Ensembl Identifiers....")
    add_MANE(processed_tables, chroms)
    print("Adding Molecular Consequences....")
    add_molecular_consequences(chroms, processed_tables)
    print("Appending HPA data.....")
    make_gene_tables(processed_tables, HPApath)
    print("Writing new HGVS Lookup table.....")
    make_HGVStable(processed_tables, chroms)


    with open(lastUpdate_file, "w") as f:
        today = date.today()
        f.write(str(today))
    f.close()


def check_updates(lastUpdate_file):
    '''
    Determines if an update is needed based on checking the date of txt file
    '''

    f = open(lastUpdate_file, "r").readlines()
    lastdate = datetime.strptime(f[0], '%Y-%m-%d').date()
    if (date.today() - lastdate).days > 31:
        updateTables()
    else:
        print('You are using the latest clinvar data')
        print(f"Last updated {lastdate}")

from datetime import datetime, date
import pandas as pd
import regex as re
import sys
from dataH import DataHandler
from Bio.Seq import Seq

#query, qtype = str(sys.argv[1]),str(sys.argv[2]) ##This needs to change to a yaml file

#test
query, qtype = 'NM_000152.5:c.271G>A','HGVS'
outfolder = '/home/thudson/edit_test_out/'

#####
# Data Query Functions
#####
def mutationType(codon1,codon2):
    aa_groups =  [["G","A","V","L","I"],
                 ["S","C","U","T","M"],
                 ["F","Y","W"],
                 ["H","K","R"],
                 ["D","E","N","Q"]]
    codon1,codon2 = Seq(codon1),Seq(codon2)
    aa1, aa2 = codon1.translate(), codon2.translate()
    mtype = ""
    if aa1== aa2:
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


def getBE(dh,guides,strand,editor,win_size= [4,8]):
    '''
    Finds codon level SNV and determines if the Base Editor Conversion can work
    *TO DO : Check HGVS term to find protein change and then bypass this 'manual' conversion
    '''

    BE_guides = {'gRNA': [],
                 'PAM' : [],
                 'SNV Position': [],
                 'Codon Position':[],
                 'Strand':[],
                 'Reference (Codon>AA)':[],
                 'Alternate (Codon>AA)': [],
                 'BE Converted (Codon>AA)': [],
                 'Conversion Type':[],
                 'Bystander': []
                 }
    hgvs_id = dh.vardf.Name.iloc[0]

    if editor == 'ABE':### A --> G
        f_delta = ['A','G']
        r_delta = ['T','C']
    else: ### C --->T
        f_delta = ['C','T']
        r_delta = ['G','A']

    #Determine Reading frame
    coding_pos = round(int(re.search("(\d+)\w+>\w+", hgvs_id).groups()[0]) /3,1) #find exon VCF position
    rf = 2 if coding_pos - round(coding_pos) == 0 else 0 if (coding_pos - round(coding_pos)) < 0.3 else 1

    for i in range(len(guides['gRNA'])):

        gseq = Seq(guides['gRNA'].iloc[i])
        SNVpos = guides.SNV_position.iloc[i]
        target_bases = Seq(gseq[win_size[0]-1:win_size[1] + 1]) #Bases inside 4-8 window
        codon_start = SNVpos - rf #which place is the SNV in the codon

        #Alternative codon
        alt_codon = gseq[codon_start - 3:codon_start] if strand == "-" else gseq[codon_start:codon_start + 3]
        aa_alt = alt_codon.translate()

        #Reference codon
        ref_codon = Seq("".join(alt_codon[x] if x!= rf else dh.ref_allele for x in [0,1,2]))
        aa_ref = ref_codon.translate()

        ##Converted allele
        convert = str(Seq(r_delta[1]) if strand == "-" else f_delta[1])
        new_codon = Seq("".join(alt_codon[x] if x != rf else convert for x in [0, 1, 2]))
        aa_new = new_codon.translate()

        mtype = mutationType(codon1 = ref_codon,codon2 = new_codon)

        ### If conversion leads to Ref change or REf change keep
        if mtype == 'Synonymous' or mtype == 'Silent' or mtype == 'Conservative':
             BE_guides['gRNA'] += [str(gseq)]
             BE_guides['PAM']+= [guides.PAM.iloc[i]]
             BE_guides['SNV Position']+= [SNVpos]
             BE_guides['Strand'] += [strand]
             BE_guides['Codon Position'] += [codon_start]
             BE_guides['Reference (Codon>AA)'] += [f"{ref_codon}>{aa_ref}"]
             BE_guides['Alternate (Codon>AA)'] += [f"{alt_codon}>{aa_alt}"]
             BE_guides['BE Converted (Codon>AA)'] += [f"{new_codon}>{aa_new}"]
             BE_guides['Conversion Type'] += [mtype]
             BE_guides['Bystander'] += [target_bases.count(f_delta[0])]

    return BE_guides

def aggGuides(dh):
    '''
    Find guides for spCas9 ABE and CBE editors

    '''
    variant_pos = dh.vardf["PositionVCF"].iloc[0]

    #find guides with NGG PAM
    guides = dh.get_Guidetable(dh.chrom, int(variant_pos) , win_size = [4,8])

    if guides.shape[0] > 0:
        guides = guides.reset_index()
        guides_dict = {}
        guides_dict['spCas9'] = guides.to_dict("list")

        cbe_guides,abe_guides = {}, {}

        if dh.alt_allele == 'C' and guides[guides.strand == '+'].shape[0] >0:
            cbe_guides = getBE(dh,guides[guides.strand == '+'], strand = '+',editor = 'CBE',win_size=[4, 8])

        if dh.alt_allele == 'G' and guides[guides.strand == '-'].shape[0] >0:
            cbe_guides = getBE(dh,guides[guides.strand == '-'],strand = '-',editor = 'CBE', win_size=[4, 8])

        if dh.alt_allele == 'A' and guides[guides.strand == '+'].shape[0] >0:
            abe_guides = getBE(dh,guides[guides.strand == '+'],strand = '+',editor = 'ABE', win_size=[4, 8])

        if dh.alt_allele == 'T' and guides[guides.strand == '-'].shape[0] >0:
            abe_guides = getBE(dh,guides[guides.strand == '-'],strand = '-',editor = 'ABE', win_size=[4, 8])

        if len(abe_guides.keys()) > 0:
            guides_dict['ABE'] = abe_guides

        if len(cbe_guides.keys()) > 0:
            guides_dict['CBE'] = cbe_guides

    else:
        guides_dict = None

    return guides_dict

#def guideCSV(dh,guides):
    #spCas9_df =
#    BE_df = pd.DataFrame(columns = ['Editor','gRNA', 'PAM','scoreDesc','fusi','crisprScan','oof'])


def summary(dh,guides, outfolder):
    '''
    text output and phenotype CSV, Guides will be outputted in a seperate CSV file
    '''
    var_out = f"{outfolder}Variant_Report_{dh.hgvs_id}.csv"
    fname = dh.prefix + "_" +re.search("[.+-]*(\d+\w+>\w+)", dh.hgvs_id).captures()[0].replace(">","_").replace(".","")
    summary_out = f"{outfolder}Summary_Report_{fname}.txt"
    #guides_out = f"{outfolder}Guides_Report_{dh.hgvs_id}.csv")
    #guideCSV(guides_out)
    vardf = dh.vardf
    dh.vardf.T.to_csv(var_out)

    variant_pos = vardf["PositionVCF"].iloc[0]
    l = len(vardf.Name.iloc[0]) + 10
    with open(summary_out, "w") as out:
        print(datetime.today().strftime("%m/%d/%y"),file = out)
        print("\n", file=out)
        print("-----", "".join("-" for x in range(l)), "-----",file = out)
        print("-----","".join(" " for x in range(l)),"-----",file = out)
        print("-----     ",vardf.Name.iloc[0],"     -----",file = out)
        print("-----", "".join(" " for x in range(l)), "-----",file = out)
        print("-----","".join("-" for x in range(l)),"-----",file = out)
        print("SNV Type: ", vardf.Type.iloc[0],file = out)
        print("Gene Name: ", vardf.GeneSymbol.iloc[0],file = out)
        print("\n", file = out)
        print("  >>>>>> Phenotype <<<<<<   ", file = out)
        print("Listed Phenotypes: ", vardf.PhenotypeList.iloc[0],file = out)
        print("Tissue Enrichment: ", vardf['RNA tissue cell type enrichment'].iloc[0],file = out)
        print("Molecular Function: ", vardf['Molecular function'].iloc[0],file = out)

        #Print Base Editors
        if len(guides.keys()) > 1:
            editor = 'ABE' if 'ABE' in guides.keys() else 'CBE'
            print("\n", file=out)
            print("  >>>>>> ABE Guides Found <<<<<< ", file=out)
            print("SNV Reference Translation: ", guides[editor]['Reference (Codon>AA)'][0], file=out)

            for i in range(len(guides[editor]['gRNA'])):
                print(f"  ---- {editor} Guide{i + 1} -----", file=out)
                print("guideSeq: ", file=out)
                print("".join(" " if p != guides[editor]['SNV Position'][i] else "*" for p in range(23)), file=out)
                print(guides[editor]['gRNA'][i], file=out)
                if guides[editor]['Strand'][i] == '+':
                    codon = "".join(" " if p != guides[editor]['Codon Position'][i] else f"ABC" for p in range(23))

                else:
                    codon = "".join(" " if p != guides[editor]['Codon Position'][i]-3 else f"ABC" for p in range(23))
                print(codon.replace("ABC", f"|{guides[editor]['Alternate (Codon>AA)'][i][-1]}|"), file=out)
                print(f"Editing Outcome   {guides[editor]['Alternate (Codon>AA)'][i][-1]} -->  {guides[editor]['BE Converted (Codon>AA)'][i][-1]}" , file=out)
                print("PAM: ", guides[editor]['PAM'][i], file=out)
                if guides[editor]['Conversion Type'][i] == 'Conservative':
                    print("** Converted codon translates is not the same amino acid as the Reference", file = out)
                print("Conversion Type: ", guides[editor]['Conversion Type'][i], file=out)
                if guides[editor]['Bystander'][i] >0:
                    print("** Warning: this editor has bystander bases in the target window", file = out)

        #Print spCas9 to file
        print("\n", file =out)
        print("  >>>>>> spCas9 Guides Found <<<<<< ", file = out)
        cnt = 0
        for i in range(len(guides['spCas9']['gRNA'])):
            print(f"  ---- spCas9 Guide{cnt+1} -----",file = out)
            print("guideSeq: ",file = out)
            print("".join(" " if p != guides['spCas9']['SNV_position'][i] else "*" for p in range(23)),file = out)
            print(guides['spCas9']['gRNA'][i],file = out)
            print("PAM: ", guides['spCas9']['PAM'][i],file = out)
            x = guides['spCas9']['Scores'][i]
            scores = re.findall("(\w+:[()%0-9]+)", x.replace(" ", ""))
            for s in scores:
                print(f"  {s}", file =out)
            cnt += 1

    out.close()
    report = open(summary_out, "r")
    for line in report:
        print(line.rstrip())



def check_HGVS(query):
    '''
    this will check the HGVS term and dertermine if it meets requirments for query
    --> check c. is present to confirm coding
    --> check if single subsitution
    ---> check not if not UTR [.c-]
    '''
    pass

def query_HGVS(hgvs_id):
    dh = DataHandler()
    vardf = dh.get_ClinVartable(hgvs_id)
    guides = aggGuides(dh)
    return dh, guides, vardf


def queryGuides(query, qtype,outfolder):
    '''
    query guides Based on the type of query 'HGVS' or 'Cell type' etc.
    ** right now limited to HGVS
    '''
    if qtype == 'HGVS':
        dh, guides, vardf = query_HGVS(hgvs_id=query)
    summary(dh,guides, outfolder)
    return guides, vardf


guides, vardf = queryGuides(query, qtype, outfolder)


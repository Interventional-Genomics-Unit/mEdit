import pandas as pd
import vcf
import os
import pickle

from dataH import DataHandler


#################################
# searches for changes in guides found in hg38 based on variants in a user submitted VCF
#################################

## assumption that the VCF imported. Has the Alternate gene specified
## A deletion greater than 5

resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"
vcf_fname = "/groups/clinical/projects/editability/tables/raw_tables/VCFs/HG02257.filtered.vcf.gz"

paths = os.listdir(resultsfolder)
hg38guide_results =[x for x in paths if x.endswith('Guides_found.csv')]
if len(hg38guide_results) > 1:
    print('warning: there but only be one file ends in Guides_Found.csv in results folder')
else:
    hg38guide_results = resultsfolder + hg38guide_results[0]
searchp_path = [resultsfolder + x for x in paths if x.endswith('guide_search_params')][0]
sitep_path = [resultsfolder + x for x in paths if x.endswith('snv_site_info')][0]
altgenome_name = 'HG02257'
refgenome_name = 'HG38'

def extract_vcf_record(snv_coord,vcf_fname,window = 30):
    #snv_coord = 'chr11:5226788'
    #extracted_seq = 'ATCCCCAAAGGACTCAAAGAACCTCTGGGTTCAAGGGTAGACCACCAGCAGCCTAAGGGT'
    ch, pos = snv_coord.split(':')
    records_found = []

    vcf_reader = vcf.Reader(filename = vcf_fname, compressed = True)
    for record in vcf_reader.fetch(ch, int(pos) - window, int(pos) + window):
        records_found.append(record)
        print(record)
    return records_found

def extract_variant_info(record,hg38extracted_seq,hg38coord):
    vt = record.var_type
    x = record.samples[0]
    alt_alleles = record.ALT
    print(record.ALT,record.REF)
    zyg = 'heterozygous' if x.is_het else 'homozygous'
    if len(alt_alleles) == 2:
        zyg = zyg + '-biallelic'
    if len(alt_alleles) > 2:
        zyg = zyg + '-multiallelic'

    alt = alt_alleles[0]
    hg38_start = int(hg38coord.split(':')[1]) - int(len(hg38extracted_seq)/2)
    vstart, vend = record.start, record.end #reference range
    rel_pos = ((vstart -hg38_start)+ 1,(vend-hg38_start)+1)
    print(rel_pos)
    print(vstart,vend)
    print(alt.type)
    if len(record.REF) == len(alt): ## substitution
        vt = vt + '-sub'
        new_seq = hg38extracted_seq[0:rel_pos[0]] + alt.sequence + hg38extracted_seq[rel_pos[1]:]
        print(new_seq)

    elif len(record.REF) > len(alt): # deletion
        vt = vt + '-del'
        #new_seq = hg38extracted_seq[0:rel_pos[0]] + alt.sequence + hg38extracted_seq[rel_pos[1]:]
        new_seq = 'undetermined'

    elif len(record.REF) < len(alt): # insertion
        vt = vt + '-ins'
        new_seq = hg38extracted_seq[0:rel_pos[0]] + alt.sequence + hg38extracted_seq[rel_pos[1]:]
        new_seq = new_seq[0:len(hg38extracted_seq)-1]
        print(new_seq)

    else:
        new_seq = 'undetermined'

    return record.REF, alt.sequence, zyg, vt, new_seq



def find_overlapping_variants(vcf_fname,sitep_path,searchp_path,altgenome_name):

    #get hg38_snv_info (60bp extracted sequence, translation info, hg38_coordinates etc.
    hg38_snvinfo = pickle.load(open(sitep_path,'rb'))

    #get guide search parameters
    search_params = pickle.load(open(searchp_path,'rb'))

    new_guides = {}
    v_info = [[],[],[],[],[],[],[]]
    for ch, data in hg38_snvinfo.items():
        for d in data:
            query, tid, eid, strand, hgref, hgalt, feature_annotation, hg38extracted_seq, codons, hg38coord = d

            #Check VCF if variant exsists in hg38 extracted_seq
            records_found = extract_vcf_record(snv_coord = hg38coord, vcf_fname = vcf_fname, window=30)

            if len(records_found) > 0:
                #Create new variant incorporated extracted seq
                for rec in records_found:
                    ref, alt, zyg, vtype, var_seq = extract_variant_info(rec,str(hg38extracted_seq),hg38coord)
                    v_info[0].append(query)
                    v_info[1].append(vtype)
                    alts = ",".join([str(x) for x in rec.samples[0].site.alleles][1:])
                    v_info[2].append(f'{ref}|{alts}')
                    v_info[3].append(alt)
                    v_info[4].append(f'chr{ch}:{rec.POS}')
                    v_info[5].append(zyg)
                    if var_seq != 'undetermined':
                        dh = DataHandler(query, strand, hgref, hgalt, feature_annotation, var_seq, codons, hg38coord)
                        new_guides_set, be_none = dh.get_Guides(search_params)
                        if len(new_guides_set['gRNA']) > 0:
                            if len(new_guides.keys()) == 0:
                                for k, v in new_guides_set.items():
                                    new_guides[k] = v
                            else:
                                for k, v in new_guides_set.items():
                                    new_guides[k] += v
                        v_info[6].append('placeholder')
                    else:
                        v_info[6].append('undetermined')
                            #print(ref, alt, zyg, vtype)
                            #print(d)
    labels = ['QueryTerm',f'{altgenome_name}_Variant_Type','REF|ALT','Examined_ALT','Var_Position', f'{altgenome_name}_Zygosity',f'{altgenome_name}_guide_impact']
    variants_found = dict(zip(labels,v_info))

    return variants_found, new_guides


def write_results(hg38guide_results,vcf_fname,sitep_path,
                  searchp_path,altgenome_name,refgenome_name):

    variants_found, new_guides_dict = find_overlapping_variants(vcf_fname,sitep_path,searchp_path,altgenome_name)

    if len(variants_found['QueryTerm']) > 1:
        hg38_gdf = pd.read_csv(hg38guide_results)
        var_df = pd.DataFrame(variants_found).set_index('QueryTerm')
        new_guides = pd.DataFrame(new_guides_dict).drop(columns=[ 'Guide_ID','SNV Position', 'Ref>Alt','Annotation'])
        hg38_gdf = hg38_gdf.loc[hg38_gdf.QueryTerm.isin(set(variants_found['QueryTerm']))]
        hg38_gdf = hg38_gdf.drop(columns =['SNV Position', 'Ref>Alt'])
        old_guides = hg38_gdf.join(var_df, how = 'outer' ,on = 'QueryTerm').reset_index(drop=True)

        new_info = new_guides.to_dict('tight')['data']
        print('n new guides', len(new_info))
        print('n old guides', len(old_guides['QueryTerm']))
        new_rows = []
        for idx, row in old_guides.iterrows():
            cnt = 1
            query, ed, coord, grna, pam = row['QueryTerm'],row['Editor'],row['Coordinates'],row['gRNA'],row['Pam']
            print( query, ed, coord, grna, pam,row['REF|ALT'])
            if row[f'{altgenome_name}_guide_impact'] == 'undetermined':  # guide never determined (ALT is deletion)
                new_rows.append(list(row[0:8]) + ['-','-','-','undetermined'] + list(row[9:14]))
                print('undetermined')
                cnt -= 1
            else:
                for n in new_info:
                    if ed == n[1]:
                        if [query,pam,grna] == [n[0],n[5],n[4]]:  # unchanged
                            print('unchanged', pam, '->', n[5], grna, '->', n[4])
                            new_info.remove(n)
                            cnt -= 1
                            break
                        elif [query,pam] == [n[0],n[5]]:  # same pam which means change is in grna
                            new_rows.append(list(row[0:8]) + [n[4], n[5], n[6], 'grna_changed_conserved'] + list(row[9:14]))
                            print('grna_changed_conserved', pam, '->', n[5], grna, '->', n[4])
                            new_info.remove(n)
                            cnt -= 1
                            break
                        elif[query,grna] == [n[0],n[4]]:
                            new_rows.append(
                                list(row[0:8]) + [n[4], n[5], n[6], 'pam_changed_conserved'] + list(row[9:14]))
                            print('pam_changed_conserved', pam, '->', n[5], grna, '->', n[4])
                            new_info.remove(n)
                            cnt -= 1
                            break
                        else:
                            pass
                if cnt < 0:  # old guides remaining and not matched == no longer exsist
                    new_rows.append(list(row[0:8]) + ['removed', 'removed', 'removed', 'pam_changed_removed'] + list(row[9:14]))
                    print('pam_changed_removed', pam, '->', '-', grna, '->', '-')

        if len(new_info) > 0:  # new guides remaining and not matched == new guides are made
            id_cnt = 1
            for n in new_info:
                x = [n[0],n[1],f'{n[1]}_NEW{id_cnt}',n[2],n[3],'-', '-', '-','-',
                                 n[4], n[5], n[6]] + list(var_df.loc[n[0]])
                new_rows.append(x[:-1])
                print('pam_changed_added', '-', '->', n[5], '-', '->', n[4])
                id_cnt += 1
        cols = ['QueryTerm', 'Editor', 'Guide_ID', 'Coordinates', 'Strand', f'{refgenome_name}_gRNA',
       f'{refgenome_name}_Pam', f'{refgenome_name}_Doench Score',f'{refgenome_name}_gRNA',
       f'{altgenome_name}_Pam', f'{altgenome_name}_Doench Score',f'{altgenome_name}_Guide_Impact',
        'Variant_Type', 'REF|ALT','Examined_ALT', f'{altgenome_name}Var_Position', f'{altgenome_name}_Zygosity']
        df = pd.DataFrame(new_rows, columns=cols)
        df.to_csv(f'{resultsfolder}/{altgenome_name}_guide_differences.csv',index = False)



'''
df = pd.read_csv("/groups/clinical/projects/editability/tables/processed_tables/variant_tables/16_variant.txt")
for x in df['HGVS_ID']:
    if x.startswith('NC'):
        print(x)
records_found = []
for p in list(df.PositionVCF):
    snv_coord = 'chr16:' + str(p)
    record = extract_vcf_record(snv_coord,vcf_fname,window = 30)
    if len(record) > 0:
        records_found.append(record)
        print(df['HGVS_ID'].loc[df['PositionVCF']==p])



'''

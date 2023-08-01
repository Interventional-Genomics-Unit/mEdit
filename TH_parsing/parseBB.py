import bbi
import subprocess
import tempfile

'''
Created 8/1/23 by TH
parses  Crispr bb and makes files per chromosome
'''

def bed_to_bigbed(df, out, btype = "bed6+3",bedToBigBed_path = "/groups/clinical/projects/editability/bedToBigBed", chrom_sizes_path ="/groups/clinical/projects/editability/hg38.chrom.sizes"):
    bed = df.copy()
    columns = bed.columns.to_list()
    bed["chrom"] = bed["chrom"].astype(str)
    bed = bed.sort_values(["chrom", "start", "end"])

    with tempfile.NamedTemporaryFile(suffix=".bed", dir = out.replace(out.split("/")[-1],"")) as f:
        bed.to_csv(
            f.name, sep="\t", columns=columns, index=False, header=False, na_rep="nan")

        cmd = f"{bedToBigBed_path} -type={btype} -tab {f.name} {chrom_sizes_path} {out}"
        p = subprocess.run(cmd,shell=True,
                           capture_output=True,
        )
        print(p.stderr)
    return p


def parse_CRISPRbb(bbi_path,n,start,end):

    with bbi.open(bbi_path) as f:

        df = f.fetch_intervals(chrom = n, start=start, end=end).drop(
            columns=["_offset", "_mouseOver", "thickStart", "thickEnd", "reserved", "_crisprScanColor",
                     "_specColor"])
        '''
        condensed terms 
        impossible to target =
            'This guide sequence is not unique in the genome. The specificity scores were not determined.'
        hard to target =
            'This guide has too many potential off-targets. The specificity score could not be calculated.'
        '''

        df['scoreDesc'] = df['scoreDesc'].replace(
            to_replace='This guide sequence is not unique in the genome. The specificity scores were not determined.',
            value='impossible to target').replace(
            to_replace='This guide has too many potential off-targets. The specificity score could not be calculated.',
            value='hard to target')
        ##combined score in dictionary
        df['ot_score'] = '{scoreDesc:' + df.scoreDesc + ' fusi:' + df.fusi + ' crisprScan:' + df.crisprScan + ' doench:' + df.doench + ' oof:' + df.oof +'}'
        df = df.drop(columns =['scoreDesc', 'fusi','crisprScan','doench','oof'])
    return df


def parse_compress_bed(bbi_path):
    chrom_size = bbi.chromsizes(bbi_path)
    chrom_names = chrom_size.keys()

    for n in chrom_names:
        df = parse_CRISPRbb(bbi_path, n = n , start = 0, end = chrom_size[n])
        out = f"/groups/clinical/projects/editability/CRISPR_hg38/{n}.bb"
        bed_to_bigbed(df=df, out=out)


bbi_path = '/groups/clinical/projects/editability/crispr.bb'
parse_compress_bed(bbi_path)

#test
#n,start,end = 'chr11_KI270721v1_random',0,100316

# Installed Modules
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio.SeqUtils import seq3
# Project Modules
import scoring
from scoring import load_model_params


class DataHandler:
    """
    search for guides given a genomic sequence and SNV info
    """

    def __init__(self, query, strand, ref, alt, feature_annotation, models_dir, extracted_seq, rf, coord,gname,dist_from_cutsite = 7):
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
        self.models_dir = models_dir
        self.coord = coord
        self.chrom = coord.split(':')[0].replace('chr', '')
        self.gname = gname
        self.dist_from_cutsite = dist_from_cutsite
        self.models_not_loaded = True
        self.models = {}

        # outputs
        self.guides_found = {'QueryTerm': [], 'GeneName':[],'Editor': [], 'Guide_ID': [], 'Coordinates': [],
                             'Strand': [], 'gRNA': [], 'Pam': [], 'Extended Guide Site': [],'Variant Position': [],
                             'Azimuth Score': [], 'DeepCas9 Score':[],'DeepCpf1 Score':[],
                             'OOF Score':[],'Ref>Alt': [], 'Annotation': []}

        self.BEguides_found = {'QueryTerm': [], 'GeneName':[],'Editor': [], 'Guide_ID': [], 'Coordinates': [],
                               'Strand': [],'gRNA': [], 'Pam': [], 'Variant Position': [],'Extended Guide Site': [],
                               'ABE score':[], 'CBE Score':[],
                               'Hg38 Reference (Codon>AA)': [], 'Alternate (Codon>AA)': [], 'BE Converted (Codon>AA)': [],
                               'Conversion Type': [], 'Bystander': [], 'Annotation': []}

    def get_score(self,site,score_name,models_dir):
        if self.models_not_loaded:
            deepcas9_model, deepcas9_sess = load_model_params('deepspcas9', models_dir)
            model1,model2 = load_model_params('deepcpf1', models_dir)
            az_model = load_model_params('azimuth', models_dir)
            self.models = {'azimuth':az_model,"deepspcas9":(deepcas9_model, deepcas9_sess),"deepcpf1":(model1,model2)}
        if score_name=='azimuth':
            site_score = scoring.azimuth(site,self.models[score_name])[0]
        if score_name == 'deepspcas9':
            model, sess_path = self.models[score_name]
            site_score = round(scoring.deepspcas9(site,model, sess_path)[0][0],2)
        if score_name == 'deepcpf1':
            model1,model2 = self.models[score_name]
            site_score = round(scoring.deepcpf1(site,model1,model2)[0][0],2)
        return site_score

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
        snv_rel_pos = len(self.extracted_seq) / 2

        for i in range(len(guides)):
            scores = {'abe': '-',
                      'cbe': '-'}
            editor, guide, pam_found, guide_strand, snvpos, editor_scores,extended_guide, start, end = guides[i]

            target_bases = Seq(guide[win_size[0] - 1:win_size[1] + 1])  # Bases inside 4-8 window

            # Converted case
            convert = str(conversion[1])
            bystander = target_bases.upper().count(conversion[0]) - 1

            if self.annotation.split(" ")[0] not in ['exon','start_codon','stop_codon']:
                ctype = 'NA'
                self.add_BEguides(name,
                                  guide,
                                  pam_found,
                                  guide_strand,
                                  snvpos,
                                  scores,
                                  extended_guide,
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
                                  scores,
                                  extended_guide,
                                  start,
                                  end,
                                  ref,
                                  alt,
                                  convert,
                                  ctype,
                                  bystander)

    def add_guides(self, name, guide, pam_found, strand, snvpos,scores,extended_guide,start,end):
        self.guides_found['QueryTerm'].append(self.query)
        self.guides_found['GeneName'].append(self.gname)
        self.guides_found['Editor'].append(name)
        self.guides_found['Guide_ID'].append(f'{name}_')
        self.guides_found['Coordinates'].append(f'{self.chrom}:{start}-{end}')
        self.guides_found['Azimuth Score'].append(scores['azimuth'])
        self.guides_found['DeepCas9 Score'].append(scores['deepcas9'])
        self.guides_found['DeepCpf1 Score'].append(scores['deepcpf1'])
        self.guides_found['OOF Score'].append(scores['oof'])
        self.guides_found['Extended Guide Site'].append(extended_guide)
        self.guides_found['Strand'].append(strand)
        self.guides_found['Pam'].append(str(pam_found))
        self.guides_found['gRNA'].append(str(guide))
        self.guides_found['Ref>Alt'].append(f"{self.NM_ref_allele}>{self.NM_alt_allele}")
        self.guides_found['Variant Position'].append(snvpos)
        self.guides_found['Annotation'].append(self.annotation)

    def add_BEguides(self, name, guide, pam_found, strand, snvpos, scores,extended_guide,start, end, ref, alt, convert, ctype, bystander):
        self.BEguides_found['QueryTerm'].append(self.query)
        self.BEguides_found['GeneName'].append(self.gname)
        self.BEguides_found['Guide_ID'].append(f'{name}_')
        self.BEguides_found['Editor'].append(name)
        self.BEguides_found['Coordinates'].append(f'{self.chrom}:{start}-{end}')
        self.BEguides_found['gRNA'].append(str(guide))
        self.BEguides_found['Pam'].append(str(pam_found))
        self.BEguides_found['Variant Position'].append(snvpos)
        self.BEguides_found['ABE score'].append(scores['abe']),
        self.BEguides_found['CBE Score'].append(scores['cbe']),
        self.BEguides_found['Strand'].append(strand)
        self.BEguides_found['Extended Guide Site'].append(extended_guide)
        self.BEguides_found['Hg38 Reference (Codon>AA)'].append(ref)
        self.BEguides_found['Alternate (Codon>AA)'].append(alt)
        self.BEguides_found['BE Converted (Codon>AA)'].append(convert)
        self.BEguides_found['Conversion Type'].append(ctype)
        self.BEguides_found['Bystander'].append(bystander)
        self.BEguides_found['Annotation'].append(self.annotation)

    def get_guide_set(self, name, pam, pamISfirst, win_size, guidelen,cut_site_position, BEmode):
        """
        :param name: editor name
        :param pam: pam seq ex:'NGG'
        :param pamISfirst: 5'or3'PAM ex:True/False
        :param win_size: list containing upper and lower limits of the targetable window. Ex:[4,8]
        :param guidelen: guide without pam length
        :return: Guide Dictionary
        """
        guides = []
        pamlen = len(pam)
        sitelen = guidelen + pamlen
        var_relative_pos = int(len(self.extracted_seq)/2)

        if BEmode:
            # Adjust window to 5' spacer
            # Example Base editing window = 4-8
            #     4   8
            #  XXXXXXXXXXXXXXXXXXXXNGGxxxxxx
            win_size = [win_size[0] - guidelen,win_size[1] -guidelen]
            pam_min, pam_max = int((var_relative_pos + abs(win_size[1]))), int((var_relative_pos + abs(win_size[0]))) + 1

        else:
            # Adjust window to 3' PAM
            # Example cas9 editing window  with cut site position -3 and dist 7= -10 and +4
            #           -10        1   4
            #  XXXXXXXXXXXXXXXXXXXXNGGxxxxxx

            x = int((var_relative_pos - abs(win_size[1])))
            pam_min, pam_max = x, (x + abs(win_size[1] - win_size[0]))
            if pamISfirst == True:
                # Adjust window to 5' PAM
                # Example cas12 editing window  with cut site position 18 and dist 7 = 15 and +29
                #  1             15             29
                #  TTTVXXXXXXXXXXXXXXXXXXXXxxxxxx
                pam_min, pam_max = pam_min - pamlen, pam_max - pamlen

        # Narrow based on guide params
        for search_strand in ["+", "-"]:
            search_seq = Seq(self.extracted_seq) if search_strand == "+" else Seq(self.extracted_seq).reverse_complement()
            snv_rel_pos = var_relative_pos
            if search_strand == "-":
                pam_start, pam_end = pam_min, pam_max - 1
                if not BEmode:
                    snv_rel_pos = var_relative_pos - 1
            else:
                pam_start, pam_end = pam_min + 1, pam_max
            try:
                pam_index = SeqUtils.nt_search(str(search_seq).upper(), pam)[1:]
            except KeyError:
                print(f"The value {pam} is not a nucleotide. Please double check your input values.")
                exit(0)

            for i in pam_index:
                if i in range(pam_start, pam_end+1):
                    scores = {'azimuth': '-',
                              'deepcas9': '-',
                              'deepcpf1': '-',
                              'oof': '-'}
                    if not pamISfirst:  # 3' PAM
                        target_start = i - guidelen
                        guide = search_seq[i - guidelen:i]
                        extended_guide = str(search_seq[target_start - 3:target_start + sitelen + 4])
                        if not BEmode:
                            snvpos = snv_rel_pos - (i + cut_site_position)
                            if snvpos >= 0:
                                snvpos += 1
                            if pam == 'NGG' and guidelen == 20:
                                #Azmith only accurate for NGG pams
                                scores['azimuth'] = self.get_score([extended_guide.upper()],'azimuth',
                                                           self.models_dir)
                                scores['deepcas9'] = self.get_score([extended_guide.upper()],'deepspcas9',
                                                           self.models_dir)
                            #oof_score use only for DSB
                            mh_score, scores['oof'] = scoring.oofscore(str(search_seq[target_start - 20:target_start + sitelen + 20]).upper())
                        else:
                            snvpos = int([guide.index(x) + 1 for x in guide if x.islower()][0])

                        pam_found = str(search_seq[i:i + pamlen])

                    else:
                        target_start = i
                        guide_start = i + pamlen
                        guide = search_seq[guide_start: guide_start + guidelen]
                        pam_found = search_seq[i:guide_start]
                        extended_guide = str(search_seq[i - 5:i + sitelen + 4])
                        snvpos = snv_rel_pos - (i + pamlen + cut_site_position)
                        if snvpos >= 0:
                            snvpos += 1
                        if 'Cas12a' in name:
                            scores['deepcpf1'] = self.get_score([extended_guide.upper()],'deepcpf1',
                                                           self.models_dir)

                    start_diff = target_start -snv_rel_pos

                    if search_strand == '-':
                        start_diff = snv_rel_pos - (target_start + sitelen)
                    start = self.SNV_chr_pos + start_diff
                    end = start + sitelen
                    
                    #print(pam_min,pam_max,name,i,snvpos, extended_guide)
                    guides.append([name, guide, pam_found, search_strand, int(snvpos), scores,extended_guide, start, end])
                    if not BEmode:
                        self.add_guides(name, guide, pam_found, search_strand, int(snvpos), scores ,extended_guide, start, end)
        return guides

    def get_Guides(self, search_params, BEsearch_params = None):
        for name, params, in search_params.items():
            try:
                pam, pamISfirst, guidelen, dsb_loc = params[0:4]
            except ValueError:
                print(f"Not enough values to unpack. Currently 'params' has the following values: {params[0:4]}")
                pam, pamISfirst, guidelen, dsb_loc = (None, None, None, None)
            try:
                win_size = [int(dsb_loc) - self.dist_from_cutsite, int(dsb_loc) + self.dist_from_cutsite]
            except TypeError:
                print(f"Cannot obtain window size based on the current attributes provided. "
                      f"Please double-check the current search parameters and try again: {params[0:4]}.")
                exit(0)

            guides = self.get_guide_set(name, pam, pamISfirst, win_size, guidelen,dsb_loc, BEmode=False)

        # if BE mode is on
        if BEsearch_params is not None:
            if len(BEsearch_params.keys()) and len(self.NC_ref_allele + self.NC_alt_allele) ==2:
                for k, params, in BEsearch_params.items():
                    pam, pamISfirst, guidelen, win_size = params[0][0], params[0][1], params[0][2], params[0][3]
                    bguides = self.get_guide_set(k, pam, pamISfirst, win_size, guidelen,cut_site_position=100, BEmode=True)

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

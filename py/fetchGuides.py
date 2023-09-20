###############
#
# Main Script with Fetch_Guides Class for running pipeline
#
###############

import regex as re
import sys
from datetime import datetime, date
import pandas as pd
from Bio.Seq import Seq
from Bio import motifs

from py.dataH import DataHandler
#from py.validate import Validator





class Fetch_Guides:


    def __init__(self, queries:list, qtype:str, editor:str, BEmode:str,resultsfolder:str,datadir:str,fasta_path:str,**kwargs):
        '''
        :param queries: list of query terms
        :param qtype: 'hgvs' or 'coord' #coord is not function yet
        :param editor: 'all', 'custom'(custom must contain kwargs, see below), or selected list or str() of the editor choices
        :param BEmode: 'off','default','all', or select BE editor for base editor choices below
        :param kwargs:

        ** if 'custom' selcted as editor in kwargs must include pam, pamISFirst, window_size (optional:name)
        ** if not custom or
        '''

        ##-----------------User Inputs--------------------##
        self.queries = queries
        self.qtype = qtype
        self.editor = editor
        self.BEmode = BEmode
        self.kwargs = kwargs

        #input paths
        self.datadir = "/groups/clinical/projects/editability/tables"
        self.processed_tables = f"{self.datadir}/processed_tables/"  # folder with cleaned clinvar/hpa tabs
        self.fasta_path = fasta_path

        #outputfolder
        self.resultsfolder = resultsfolder

        ##---------------libraries and keys--------------------##
        # [editor]: pam, pamISfirst, win_size, guidelen, scoring, notes/altnames
        self.editor_choices = ['spCas9', 'saCas9', 'spG', 'SpRY-HighE', 'LbCpf1',
                               'scCas9', 'stCas9', 'iSpyMacCas9', 'CasX', 'Cpf1']

        self.BE_choices = ['BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-Gam', 'YE1-BE3', 'EE-BE3', 'YE2-BE3',
                           'YEE-BE3','VQR-BE3','VRER-BE3','SaBE3','SaBE4','SaBE4-Gam','Sa(KKH)-BE3','xBE3','eA3A-BE3',
                           'A3A-BE3','BE-PLUS','ABE7.9','ABE7.10','xABE,ABESa','VQR-ABE','VRER-ABE','Sa(KKH)-ABE']



        self.editor_pamlib = {'spCas9':('NGG', False,20,[4,8], 'Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
                              'saCas9': ('NNGRRT', False, 21,[4,8], 'Cas9 S. Aureus 21 base guide'),
                              'spG':('NGN',False,20,[4,8],'20bp-NGN - SpG'),
                              'SpRY-HighE': ('NRN',False,20,[4,8],'High Efficiency Pam'),
                              'LbCpf1': ('TTTA', True, 23,[4,8], 'LbCpf1'),
                              'scCas9':('NNGT',False,20,[4,8],'20bp-NNGT - Cas9 S. canis - high efficiency PAM, recommended'),
                              'stCas9':('NNAGAA',False,20,[4,8],'Cas9 S. Thermophilus'),
                              'iSpyMacCas9':('NAA',False,20,[4,8],''),
                              'CasX':('TTCN',True,20,[4,8],'Cas12e'),
                              'Cpf1':('TTTV',True,23,[4,8],'TTT(A/C/G)-23bp - Cas12a (Cpf1)')
                              }


        ## ------------Defaults and settings------------------##
        ##configure editor options
        self.search_params = self.configure_search_params()

        #configure BE options
        if self.BEmode != 'off':
            self.BE_search_params = self.set_BE_params()


        #---------------Ouputs--------------------------##
        self.all_variant = pd.DataFrame()
        self.all_gene = pd.DataFrame()
        self.all_guides = {}
        self.all_BE = {}

        self.found_genes = []

    def set_guidelen(self,guidelen):
        self.guidelen = guidelen

    def set_win_size(self,win_size):
        self.win_size = win_size

    def configure_search_params(self):
        '''
        set paramteres for the selected editor or editors(not BE editors)
        '''
        # search for all guides
        if 'all' == self.editor:
            self.search_params = self.editor_pamlib

        #search for selected subset
        if type(self.editor) is list:
            self.search_params = {}
            for e in self.editor:
                self.search_params[e] = self.editor_pamlib[e]

        # else use single set parameters
        else:
            # default - spCas9 params
            self.win_size = [4, 8]
            self.guidelen = 20
            self.scoring = 'doench'
            self.pamISfirst = False
            self.pam = 'NGG'
            self.name = 'spCas9'
            self.notes = 'none'

            # set custom editor params
            if 'custom' == self.editor:
                self.search_params = self.set_params(self.kwargs)

            # select a single editor
            if self.editor in self.editor_choices:
                self.search_params = self.set_params({'name': self.editor,
                                                      'pam': self.editor_pamlib[self.editor][0],
                                                      'pamISfirt': self.editor_pamlib[self.editor][1],
                                                      'win_size': self.editor_pamlib[self.editor][2],
                                                      'guidelen': self.editor_pamlib[self.editor][3],
                                                      'scoring': self.editor_pamlib[self.editor][4],
                                                      'notes': self.editor_pamlib[self.editor][5]})
        return self.search_params


    def set_params(self,kwargs):
        #opts = ['pam', 'pamISfirst', 'guidelen','win_size', 'name','scoring']

        if 'pam' in kwargs.keys():
            self.pam = kwargs['pam']
        if 'pamISfirst' in kwargs.keys():
            self.pamISfirst = kwargs['pamISfirst']
        if 'guidelen' in kwargs.keys():
            self.guidelen = kwargs['guidelen']
        if 'name' in kwargs.keys():
            self.name = kwargs['name']
        if 'win_size' in kwargs.keys():
            self.win_size = kwargs['win_size']
        if 'notes' in kwargs.keys():
            self.notes = kwargs['notes']

        params = {self.name: (self.pam, self.pamISfirst, self.win_size, self.guidelen)}

        return params

    def set_multi_params(self):
        #TODO: make automatic function using configureMultipams
        pass


    def set_BE_params(self):
        # sets base editor search params, each key is a list of 2 or more; refernce seq search params,
        # then any set that follows starts with the conversion (ex. 'AG' is A --> G) and then the base editors that have the same params

        BE_lib = {'spCas9-def': [('NGG', False, 20,[4,8]),('CT','BE3', 'BE4', 'BE4max', 'BE4-Gam'),('AG','ABE7.9','ABE7.10','ABEmax')],
                       'spCas9-YE': [('NGG', False, 20, [5, 6]),('CT','YE1-BE3','EE-BE3', 'YE2-BE3','YEE-BE3')],
                       'spCas9-VQR': [('NGA', False, 20, [4, 8]), ('CT','VQR-BE3'),('AG','VQR-ABE')],
                       'spCas9-VRER': [('NGCG', False, 20, [4, 8]), ('CT','VRER-BE3'),('AG','VRER-ABE')],
                       'saCas9-KKH': [('NNGRRT', False, 20, [3, 12]),('CT','SaBE3','SaBE4','SaBE4-Gam,Sa(KKH)-BE3'),('AG','ABESa','Sa(KKH)-ABE')]}

        if self.BEmode == 'default':
            self.BE_search_params = {'spCas9-def':BE_lib['spCas9-def']}

        if self.BEmode =='all':
            self.BE_search_params = BE_lib

        if self.BEmode in self.BE_choices:
            for k,v in BE_lib.values():
                if self.BEmode in v[1]:
                     self.BE_search_params = [BE_lib[k]('CT',self.BEmode)('AT')]
                if self.BEmode in v[2]:
                    self.BE_search_params = [BE_lib[k]('CT')('AT',self.BEmode)]

        return self.BE_search_params

    def write_guide_csv(self,guides,gtype):
        df = pd.DataFrame(guides)
        nameout = 'BaseEditors_found.csv' if gtype == 'BE' else 'Guides_found.csv'
        datenow = date.today().strftime('%Y-%m-%d')
        out = f'{self.resultsfolder}_{datenow}_{nameout}'
        df.to_csv(out,index = False)
        return df

    def add_clininfo(self,clininfo):
        variantinfo = clininfo[
            ['HGVS_ID', 'GeneSymbol', 'Chr', 'Start', 'End', 'Coding Strand', 'AltAlleleVCF', 'RefAlleleVCF',
             'AlleleID', 'Type', 'GeneID', 'HGNC_ID', 'ClinicalSign',
             'ClinSigSimple', 'nsv/esv (dbVar)', 'RCVaccession', 'PhenoList', 'OriginSimple', 'ChrAccession',
             'VariationID', 'PositionVCF', 'OMIM']]
        geneinfo = clininfo[['GeneSymbol', 'Chr', 'Start', 'End', 'Coding Strand', 'Biological process',
                             'Disease involvement', 'GeneID', 'OMIM', 'IDs', 'Gene_Type', 'Molecular function',
                             'Ensembl', 'Gene synonym',
                             'Gene description', 'RNA tissue specific nTPM', 'RNA tissue distribution',
                             'RNA tissue cell type enrichment', 'RNA single cell type specific nTPM',
                             'RNA tissue specific nTPM.1']]

        if len(self.all_variant.index) == 0:
            self.all_variant = variantinfo
            self.all_gene = geneinfo

        else:
            self.all_variant = pd.concat([self.all_variant, variantinfo])

            if clininfo['GeneSymbol'].iloc[0] not in self.found_genes:
                self.all_gene = pd.concat([self.all_gene, geneinfo])
                self.found_genes.append(geneinfo['GeneSymbol'].iloc[0])

    def run_FetchGuides(self):

        for query in self.queries:
            #dh = DataHandler(fg.datadir,fg.fasta_path)
            #clininfo = dh.get_ClinVartable(queries[0], fg.qtype)
            dh = DataHandler(self.datadir,self.fasta_path)
            clininfo = dh.get_ClinVartable(query, self.qtype)

            try: # If query found in database search guides
                variant_pos = dh.vardf["PositionVCF"].iloc[0]
                if self.BEmode != 'off':
                    #guides, BEguides = dh.get_Guides(fg.search_params, fg.BE_search_params)
                    guides, BEguides = dh.get_Guides(self.search_params, self.BE_search_params)

                    if len(BEguides['gRNA']) > 0:
                        if len(self.all_BE.keys()) == 0:
                            for k, v in BEguides.items():
                                self.all_BE[k] = v
                        else:

                            for k, v in BEguides.items():
                                self.all_BE[k] += v

                else:
                    guides = dh.get_Guides(self.search_params)

                if len(guides['gRNA']) > 0:
                    if len(self.all_guides.keys()) == 0:
                        for k, v in guides.items():
                            self.all_guides[k] = v

                    else:

                        for k, v in guides.items():
                            self.all_guides[k] += v

                    print(len((guides['gRNA'])), ' guides found for ', query)
                    self.add_clininfo(clininfo)
                else:
                    print(f"No guides found for the query {query}")
            except:
                pass
        guidedf = self.write_guide_csv(self.all_guides,gtype='guides')
        BEdf = self.write_guide_csv(self.all_BE, gtype='BE')
        return self.all_variant,self.all_gene, guidedf, BEdf




#Test

#Paths---------------------------
input_file = "/groups/clinical/projects/editability/medit_queries/medit_test/test_in/hgvs_test_queries.csv"
resultsfolder = "/groups/clinical/projects/editability/medit_queries/medit_test/test_out/"
datadir = "/groups/clinical/projects/editability/tables/"
fasta_path ="/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa"

#Input Extraction-------------------
df = pd.read_csv(input_file)
queries = list(df.iloc[:,0])
qtype = 'hgvs'
BEmode = 'default'
editor = 'all'


# Get query items
fg = Fetch_Guides(queries,qtype,editor,BEmode,resultsfolder,datadir,fasta_path)
all_clin_info, all_gene, all_guides, all_BE = fg.run_FetchGuides()





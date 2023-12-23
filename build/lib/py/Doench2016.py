# from sklearn.ensemble import GradientBoostingRegressor
import pickle
import numpy as np
import itertools
import Bio.SeqUtils.MeltingTemp as Tm
import pandas as pd

### Python 2.7 version

#Fusi/doench/azmith - ontarget score, he efficiency score tries to predict --> is determined by off-targets
    #   if a guide leads to rather strong or weak cleavage. According to (Haeussler et al. 2016),
    #  the Doench 2016 Efficiency score should be used to select the guide with the highest cleavage
    # efficiency when expressing guides from RNA PolIII Promoters such as U6. Scores are given as percentiles,
    # e.g. "70%" means that 70% of mammalian guides have a score equal or lower than this guide. The raw score number is also shown in parentheses after the percentile.
#off- target scores

file = "pkl/Doench_2016_18.01_model_nopos.pickle"
f = open(file, 'rb')
model = pickle.load(f) #, encoding='bytes')
model, learn_options = model
'''
model = GradientBoostingRegressor(alpha=0.5, init=None, learning_rate=0.1, loss='ls',
             max_depth=3, max_features=None, max_leaf_nodes=None,
             min_samples_leaf=1, min_samples_split=2,
             min_weight_fraction_leaf=0.0, n_estimators=100,
             random_state=1, subsample=1.0, verbose=0,
             warm_start=False) #presort='auto
learn_options = {'include_strand': False, 'weighted': None, 'algorithm_hyperparam_search': 'grid', 'num_thread_per_proc':
             None, 'extra pairs': False, 'gc_features': True,
                 'test_genes': np.array([u'CD5', u'CD45', u'THY1', u'H2-K', u'CD28', u'CD43', 'CD33', 'CD13',
       'CD15', u'HPRT1', u'CCDC101', u'MED12', u'TADA2B', u'TADA1',
       u'CUL3', u'NF1', u'NF2'], dtype=object),
                 'testing_non_binary_target_name': 'ranks',
                 'train_genes': np.array([u'CD5', u'CD45', u'THY1', u'H2-K', u'CD28', u'CD43', 'CD33', 'CD13',
       'CD15', u'HPRT1', u'CCDC101', u'MED12', u'TADA2B', u'TADA1',
       u'CUL3', u'NF1', u'NF2'], dtype=object),
                 'cv': 'gene',
                 'adaboost_alpha': 0.5,
                 'all pairs': False,
                 'binary target name': 'score_drug_gene_threshold',
                 'normalize_features': False,
                 'nuc_features': True, 'include_gene_effect': False,
                 'num_genes_remove_train': None, 'include_gene_guide_feature': 0,
                 'include_known_pairs': False, 'include_gene_feature': False,
                 'training_metric': 'spearmanr', 'num_proc': 4, 'include_drug': False,
                 'include_microhomology': False, 'V': 3, 'include_Tm': True, 'adaboost_loss': 'ls', 'rank-transformed target name': 'score_drug_gene_rank', 'include_pi_nuc_feat': True, 'include_sgRNAscore': False, 'adaboost_CV': False, 'flipV1target': False, 'include_NGGX_interaction': True, 'seed': 1,
                 'NDGC_k': 10, 'raw target name': None, 'all_genes': np.array([u'CD5', u'CD45', u'THY1', u'H2-K', u'CD28', u'CD43', 'CD33', 'CD13',
       'CD15', u'HPRT1', u'CCDC101', u'MED12', u'TADA2B', u'TADA1',
       u'CUL3', u'NF1', u'NF2'], dtype=object), 'order': 2, 'include_gene_position': False}
'''


def countGC(s):
    return len(s[5:25].replace('A', '').replace('T', ''))


def nucleotide_features(s, order, include_pos_independent, max_index_to_use, prefix="", feature_type='all'):
    '''
    compute position-specific order-mer features for the 4-letter alphabet
    (e.g. for a sequence of length 30, there are 30*4 single nucleotide features
          and (30-1)*4^2=464 double nucleotide features
    '''
    if max_index_to_use is not None:
        s = s[:max_index_to_use]
    # assert(len(s)==30, "length not 30")
    # s = s[:30] #cut-off at thirty to clean up extra data that they accidentally left in, and were instructed to ignore in this way
    raw_alphabet = ['A', 'T', 'C', 'G']
    alphabet = ["".join(i) for i in itertools.product(raw_alphabet, repeat=order)]
    features_pos_dependent = np.zeros(len(alphabet) * (len(s) - (order - 1)))
    features_pos_independent = np.zeros(np.power(len(raw_alphabet), order))

    # for position in range(0, len(s)-order, 1): JENN 9/4/2014 failing when len(s)=2
    for position in range(0, len(s) - order + 1, 1):
        nucl = s[position:position + order]
        features_pos_dependent[alphabet.index(nucl) + (position * len(alphabet))] = 1.0
        features_pos_independent[alphabet.index(nucl)] += 1.0
    index_dependent = ['%s_pd.Order%d_P%d' % (prefix, order, i) for i in range(len(features_pos_dependent))]

    if feature_type == 'all' or feature_type == 'pos_independent':
        if include_pos_independent:
            index_independent = ['%s_pi.Order%d_P%d' % (prefix, order, i) for i in range(len(features_pos_independent))]
            if feature_type == 'all':
                return pd.Series(features_pos_dependent, index=index_dependent), pandas.Series(
                    features_pos_independent, index=index_independent)
            else:
                return pd.Series(features_pos_independent, index=index_independent)

    if np.any(np.isnan(features_pos_dependent)): raise Exception("found nan features in features_pos_dependent")
    if np.any(np.isnan(features_pos_independent)): raise Exception("found nan features in features_pos_independent")

    return pd.Series(features_pos_dependent, index=index_dependent)


def Tm_feature(data):
    '''
    assuming '30-mer'is a key
    get melting temperature features from:
        0-the 30-mer ("global Tm")
        1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA, i.e. the 5nts immediately proximal of the NGG PAM
        2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
        3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
    '''
    sequence = data['30mer'].values
    featarray = np.ones((sequence.shape[0],4))
    for i, seq in enumerate(sequence):
        if seq[25:27]!="GG":
            raise Exception("expected GG but found %s" % seq[25:27])
        rna = False
        featarray[i,0] = Tm.Tm_NN(seq)        #30mer Tm
        featarray[i,1] = Tm.Tm_NN(seq[20:25]) #5nts immediately proximal of the NGG PAM
        featarray[i,2] = Tm.Tm_NN(seq[12:20])   #8-mer
        featarray[i,3] = Tm.Tm_NN(seq[7:12])      #5-mer

    feat = pd.DataFrame(featarray, index=data.index, columns=["Tm global_%s" % rna, "5mer_end_%s" %rna, "8mer_middle_%s" %rna, "5mer_start_%s" %rna])

    return feat


def NGGX_interaction_feature(data):
    '''
    assuming 30-mer, grab the NGGX _ _ positions, and make a one-hot
    encoding of the NX nucleotides yielding 4x4=16 features
    '''
    sequence = data['30mer'].values
    feat_NX = pd.DataFrame()
    # check that GG is where we think
    for seq in sequence:
        if seq[25:27] != "GG":
            raise Exception("expected GG but found %s" % seq[25:27])
        NX = seq[24]+seq[27]
        NX_onehot = nucleotide_features(NX,order=2, include_pos_independent=False, max_index_to_use=2, prefix="NGGX")
        # NX_onehot[:] = np.random.rand(NX_onehot.shape[0]) ##TESTING RANDOM FEATURE
        feat_NX = pd.concat([feat_NX, NX_onehot], axis=1)
    return feat_NX.T


def check_feature_set_dimensions(feature_sets):
    '''
    Ensure the # of people is the same in each feature set
    '''
    N = None
    for ft in list(feature_sets.keys()):
        N2 = feature_sets[ft].shape[0]
        if N is None:
            N = N2
        else:
            assert N == N2, "# of individuals do not match up across feature sets"


def r2_DoenchScore(input_seqs,model,learn_options):
    res = []
    for seq in input_seqs:
        if "N" in seq:
            res.append(-1)  # can't do Ns
            continue

        pam = seq[25:27]
        if pam != "GG":
            # res.append(-1)
            # continue
            seq = list(seq)
            seq[25] = "G"
            seq[26] = "G"
            seq = "".join(seq)
        res.append(seq)
    seqs = np.array(res)

    learn_options["V"] = 2

    Xdf = pd.DataFrame(columns=['30mer', 'Strand'],
                        data=zip(seqs, np.repeat('NA', seqs.shape[0])))
    gene_position = pd.DataFrame(columns=['Percent Peptide', 'Amino Acid Cut position'],
                                             data=zip(np.ones(seqs.shape[0]) * -1,
                                                      np.ones(seqs.shape[0]) * -1))
    #Xdf = pd.DataFrame(columns=['30mer', 'Strand'], data=[[seqs, 'NA']])
    #gene_position = (pd.DataFrame(columns=['Percent Peptide', 'Amino Acid Cut position'],
    #                          data=[[percent_peptide, aa_cut]]))
    all_lens = Xdf['30mer'].apply(len).values
    feature_sets = {}
    pam_audit = True
    length_audit = False

    if learn_options["nuc_features"]:
    # spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
        seq_data_frame = Xdf['30mer']
        maxorder = learn_options["order"]
        max_index_to_use = 30
        for order in range(1, maxorder + 1):
            # print "\t\tconstructing order %s features" % order
            num_proc = learn_options["num_proc"]
            include_pos_independent = True,
            max_index_to_use = max_index_to_use
            prefix = ""
            feat_pd = seq_data_frame.apply(nucleotide_features, args=(
            order, include_pos_independent, max_index_to_use, prefix, 'pos_dependent'))
            feat_pi = seq_data_frame.apply(nucleotide_features, args=(
            order, include_pos_independent, max_index_to_use, prefix, 'pos_independent'))

            feature_sets['%s_nuc_pd_Order%i' % (prefix, order)] = feat_pd
            if learn_options['include_pi_nuc_feat']:
                feature_sets['%s_nuc_pi_Order%i' % (prefix, order)] = feat_pi

    if learn_options["gc_features"]:
        gc_count = Xdf['30mer'].apply(countGC)
        gc_count.name = 'GC count'
        gc_above_10 = (gc_count > 10) * 1
        gc_above_10.name = 'GC > 10'
        gc_below_10 = (gc_count < 10) * 1
        gc_below_10.name = 'GC < 10'
        feature_sets['gc_above_10'] = pd.DataFrame(gc_above_10)
        feature_sets['gc_below_10'] = pd.DataFrame(gc_below_10)
        feature_sets['gc_count'] = pd.DataFrame(gc_count)

    if learn_options["include_NGGX_interaction"]:
        feature_sets["NGGX"] = NGGX_interaction_feature(Xdf)

    if learn_options["include_Tm"]:
        feature_sets["Tm"] = Tm_feature(Xdf)

    check_feature_set_dimensions(feature_sets)

    F = feature_sets[list(feature_sets.keys())[0]].shape[0]
    for set in list(feature_sets.keys()):
        F2 = feature_sets[set].shape[0]
        assert F == F2, "not same # individuals for features %s and %s" % (list(feature_sets.keys())[0], set)

    N = feature_sets[list(feature_sets.keys())[0]].shape[0]
    inputs = np.zeros((N, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for set in list(feature_sets.keys()):
        inputs_set = feature_sets[set].values
        dim[set] = inputs_set.shape[1]
        dimsum += dim[set]
        inputs = np.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[set].columns.tolist())
    scores = model.predict(inputs)
    return scores


# Test
# cas9 guides ADA chr20:44,621,066-44,621,088
exgrna6, grnapam6 = 'TCATCTGGTAATCAGTGTCCAGGGTGGACT', 'CTGGTAATCAGTGTCCCGGGTGG'
#chr20:44621059-44621088
exgrna4, grnapam4 = 'TGGTCATCTGGTAATCAGTGTCCAGGGTGG','CATCTGGTAATCAGTGTCCC'
seqs = [exgrna4,exgrna6]
#[60,71]
scores = r2_DoenchScore(seqs,model,learn_options)
print(scores)

import math
from math import exp
import regex as re
import os
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
import pickle
import numpy as np
from featurization import featurize_data


def load_model_params(score_name):
    try:
        if score_name == 'cfd':
            mm_scores = pickle.load(open(os.path.dirname(os.path.realpath(__file__)) + '/pkl/mismatch_score.pkl', 'rb'))
            pam_scores = pickle.load(open(os.path.dirname(
                os.path.realpath(__file__)) + '/pkl/PAM_scores.pkl', 'rb'))
            return mm_scores, pam_scores
        if score_name == 'azimuth':
            model = pickle.load(open(os.path.dirname(
                os.path.realpath(__file__)) + '/pkl/python3_V3_model_no.pos.pickle', 'rb'))
            return model
        if score_name == 'doench14':
            doench14params = pickle.load(open(os.path.dirname(
                os.path.realpath(__file__)) + '/pkl/doench14params.pkl', 'rb'))
        return doench14params
    except:
        raise Exception(f"pickle file containing {score_name} could not be opened")

def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

def doench2014(cas9_sites):
        """
        Doench 2014 'on-target score'
        Input is a 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'
        Doench et al, Nat Biotech 2014, PMID 25184501, http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
        Code from crispor - https://github.com/maximilianh/crisporWebsite
        """
        doench2014params = load_model_params(score_name = 'doench14')
        scores = []
        global gcWeight
        intercept = 0.59763615
        gcHigh = -0.1665878
        gcLow = -0.2026259

        for seq in cas9_sites:
            assert (len(seq) == 30)
            score = intercept
            guideSeq = seq[4:24]
            gcCount = guideSeq.count("G") + guideSeq.count("C")

            if gcCount <= 10:
                gcWeight = gcLow
            if gcCount > 10:
                gcWeight = gcHigh

            score += abs(10 - gcCount) * gcWeight

            for pos, modelSeq, weight in doench2014params:
                subSeq = seq[pos:pos + len(modelSeq)]
                if subSeq == modelSeq:
                    score += weight
            scores.append(int(100 * (1.0 / (1.0 + math.exp(-score)))))

        return score

def azimuth(cas9_sites):
    '''
    Doench/Fusi 2016 Rule -2 on-target / efficiency score now packaged as 'Azimuth'
    This script is copied and modified to suit from https://github.com/MicrosoftResearch/Azimuth
    predicts whether a guide exhibits strong or weak cleavage
    Score range 0-100. A score higher than 55% is recommended
    '''
    #Doench/Fusi 2016 Rule -2 'on-target score'
    # This script is copied and modified to suit from https://github.com/MicrosoftResearch/Azimuth
    model = load_model_params(score_name='azimuth')
    model, learn_options = model

    res = []
    for seq in cas9_sites:
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

    feature_sets = featurize_data(data=Xdf, learn_options=learn_options, Y=pd.DataFrame())
    keys = list(feature_sets.keys())
    N = feature_sets[list(keys)[0]].shape[0]
    inputs = np.zeros((N, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for set in keys:
        inputs_set = feature_sets[set].values
        dim[set] = inputs_set.shape[1]
        dimsum += dim[set]
        inputs = np.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[set].columns.tolist())
    scores = model.predict(inputs)
    scores = [(s * 100).round(2) for s in scores]
    return scores


def oofscore(seq):
    '''
    copied and adapted code from Bae et al. https://www.nature.com/articles/nmeth.3015
    computes both microhomology and out-of-frame score
    A measurement of how likely an out-of-frame deletion occurs after a knock-out experiment
    based on microhomology
    scoring range 0-100. The higher the oof score, the more deletions have a length that is not a multiple of three
     A score above 66 is recommended
    The higher the oof score, the more deletions have a length that is not a multiple of three
    '''
    length_weight = 20.0
    left = 30
    right = len(seq) - int(left)

    s1 = []
    for k in range(2, left)[::-1]:
        for j in range(left, left + right - k + 1):
            for i in range(0, left - k + 1):
                if seq[i:i + k] == seq[j:j + k]:
                    length = j - i
                    s1.append(seq[i:i + k] + '\t' + str(i) + '\t' + str(i + k) + '\t' + str(j) + '\t' + str(j + k) + '\t' + str(length))

    if s1 != "":
        list_f1 = s1
        sum_score_3 = 0
        sum_score_not_3 = 0

        for i in range(len(list_f1)):
            n = 0
            score_3 = 0
            score_not_3 = 0
            line = list_f1[i].split('\t')
            scrap = line[0]
            left_start = int(line[1])
            left_end = int(line[2])
            right_start = int(line[3])
            right_end = int(line[4])
            length = int(line[5])

            for j in range(i):
                line_ref = list_f1[j].split('\t')
                left_start_ref = int(line_ref[1])
                left_end_ref = int(line_ref[2])
                right_start_ref = int(line_ref[3])
                right_end_ref = int(line_ref[4])

                if (left_start >= left_start_ref) and (left_end <= left_end_ref) and (
                        right_start >= right_start_ref) and (right_end <= right_end_ref):
                    if (left_start - left_start_ref) == (right_start - right_start_ref) and (
                            left_end - left_end_ref) == (right_end - right_end_ref):
                        n += 1
                else:
                    pass

            if n == 0:
                if (length % 3) == 0:
                    length_factor = round(1 / exp((length) / (length_weight)), 3)
                    num_GC = len(re.findall('G', scrap)) + len(re.findall('C', scrap))
                    score_3 = 100 * length_factor * ((len(scrap) - num_GC) + (num_GC * 2))

                elif (length % 3) != 0:
                    length_factor = round(1 / exp((length) / (length_weight)), 3)
                    num_GC = len(re.findall('G', scrap)) + len(re.findall('C', scrap))
                    score_not_3 = 100 * length_factor * ((len(scrap) - num_GC) + (num_GC * 2))
            sum_score_3 += score_3
            sum_score_not_3 += score_not_3

        mh_score = round(sum_score_3 + sum_score_not_3,2)
        oof_score = round((sum_score_not_3) * 100 / (sum_score_3 + sum_score_not_3),2)

    return mh_score, oof_score


def cfd_score(seq1, seq2,pam,mm_scores,pam_scores):
    '''
    Doench 2016 off-target scoring
    Doench, Fusi, et al.  Nature Biotechnology 34, 184–191 (2016)."

    '''
    score = 1
    shorter, longer = sorted([seq1, seq2], key=len)
    for i in range(-len(shorter), 0):
        if (seq1[i] != seq2[i]):
            key = 'r' + seq1[i] + ':d' + revcom(seq2[i]) + ',' + str(20 + i + 1)
            score *= mm_scores[key]
    score *= pam_scores[pam[-2:]]
    return score

def cfd_multi_score(spacer, off_target_seqs):
    '''
    on_target_seq is spacer site and off_target includes pam
    '''
    mm_scores, pam_scores = load_model_params(score_name='cfd')
    scores = []
    seq1 = spacer.upper().replace('T', 'U')
    sum_scores = 0

    for off_target_seq in off_target_seqs:
        m_off = re.search('[^ATCG]',off_target_seq)
        if (m_off is None):
            pam = off_target_seq[-2:]
            seq2 = off_target_seq[:20].upper().replace('T','U')
            if off_target_seq != spacer:
                score = cfd_score(seq1, seq2, pam,mm_scores,pam_scores)
                sum_scores += score
            else:
                score = -1
        scores.append(round(score,2))
        guide_cfd_score = 100 / (100+sum_scores)
        guide_cfd_score = int(round(guide_cfd_score*100))
    return scores, guide_cfd_score



'''
Test
pam = 'AGG'
spacer = 'GTGCGGCTGGCCCAGGACCT'
target = spacer + pam
target30mer = 'CTTGTGCGGCTGGCCCAGGACCTAGGCGAG'
target60mer = 'ATCTCTTACAACGACTTCTTGTGCGGCTGGCCCAGGACCTAGGCGAGGCAGTAGGGGATGACA'
off_target_seqs = ['GTGGGGCTGACCCAGGACCTGAG','GGGCCTCTGGCCCAGGACCTGGG']


print('Microhomology Out-Of-Frame score: ',oofscore(seq = target60mer))
# Ans: 3730.1, 64.2

print('Azimuth on-target score: ',azimuth([target30mer])[0])
#Ans: 0.38

print('CFD scores, Guide specifity score: ',cfd_multi_score(spacer=spacer, off_target_seqs = off_target_seqs))
#Ans: [0.12, 0.44], 99
'''


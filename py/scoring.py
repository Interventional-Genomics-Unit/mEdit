import math
from math import exp
from re import findall
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
import pickle
import numpy as np
from featurization import featurize_data

### inputs
model_file = "/home/thudson/projects/editability/pkl/python3_V3_model_no.pos.pickle"
model = pickle.load(open(model_file, "rb"))


doench14params = [
    # Doench: Doench et al, Nat Biotech 2014, PMID 25184501, http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
    # pasted/typed table from PDF and converted to zero-based positions
    (1, 'G', -0.2753771), (2, 'A', -0.3238875), (2, 'C', 0.17212887), (3, 'C', -0.1006662),
    (4, 'C', -0.2018029), (4, 'G', 0.24595663), (5, 'A', 0.03644004), (5, 'C', 0.09837684),
    (6, 'C', -0.7411813), (6, 'G', -0.3932644), (11, 'A', -0.466099), (14, 'A', 0.08537695),
    (14, 'C', -0.013814), (15, 'A', 0.27262051), (15, 'C', -0.1190226), (15, 'T', -0.2859442),
    (16, 'A', 0.09745459), (16, 'G', -0.1755462), (17, 'C', -0.3457955), (17, 'G', -0.6780964),
    (18, 'A', 0.22508903), (18, 'C', -0.5077941), (19, 'G', -0.4173736), (19, 'T', -0.054307),
    (20, 'G', 0.37989937), (20, 'T', -0.0907126), (21, 'C', 0.05782332), (21, 'T', -0.5305673),
    (22, 'T', -0.8770074), (23, 'C', -0.8762358), (23, 'G', 0.27891626), (23, 'T', -0.4031022),
    (24, 'A', -0.0773007), (24, 'C', 0.28793562), (24, 'T', -0.2216372), (27, 'G', -0.6890167),
    (27, 'T', 0.11787758), (28, 'C', -0.1604453), (29, 'G', 0.38634258), (1, 'GT', -0.6257787),
    (4, 'GC', 0.30004332), (5, 'AA', -0.8348362), (5, 'TA', 0.76062777), (6, 'GG', -0.4908167),
    (11, 'GG', -1.5169074), (11, 'TA', 0.7092612), (11, 'TC', 0.49629861), (11, 'TT', -0.5868739),
    (12, 'GG', -0.3345637), (13, 'GA', 0.76384993), (13, 'GC', -0.5370252), (16, 'TG', -0.7981461),
    (18, 'GG', -0.6668087), (18, 'TC', 0.35318325), (19, 'CC', 0.74807209), (19, 'TG', -0.3672668),
    (20, 'AC', 0.56820913), (20, 'CG', 0.32907207), (20, 'GA', -0.8364568), (20, 'GG', -0.7822076),
    (21, 'TC', -1.029693), (22, 'CG', 0.85619782), (22, 'CT', -0.4632077), (23, 'AA', -0.5794924),
    (23, 'AG', 0.64907554), (24, 'AG', -0.0773007), (24, 'CG', 0.28793562), (24, 'TG', -0.2216372),
    (26, 'GT', 0.11787758), (28, 'GG', -0.69774)]

def doench2014(cas9_sites,doench2014params):
        """
        Doench 2014 'on-target score'
        Input is a 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'
        Doench et al, Nat Biotech 2014, PMID 25184501, http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
        Code from crispor - https://github.com/maximilianh/crisporWebsite
        """
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

def azimuth(cas9_sites,model = model):
    '''
    Doench/Fusi 2016 Rule -2 on-target / efficiency score now packaged as 'Azimuth'
    This script is copied and modified to suit from https://github.com/MicrosoftResearch/Azimuth
    predicts whether a guide exhibits strong or weak cleavage
    Score range 0-100. A score higher than 55% is recommended
    '''
    #Doench/Fusi 2016 Rule -2 'on-target score'
    # This script is copied and modified to suit from https://github.com/MicrosoftResearch/Azimuth
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
                    num_GC = len(findall('G', scrap)) + len(findall('C', scrap))
                    score_3 = 100 * length_factor * ((len(scrap) - num_GC) + (num_GC * 2))

                elif (length % 3) != 0:
                    length_factor = round(1 / exp((length) / (length_weight)), 3)
                    num_GC = len(findall('G', scrap)) + len(findall('C', scrap))
                    score_not_3 = 100 * length_factor * ((len(scrap) - num_GC) + (num_GC * 2))
            sum_score_3 += score_3
            sum_score_not_3 += score_not_3

        mh_score = round(sum_score_3 + sum_score_not_3,2)
        oof_score = round((sum_score_not_3) * 100 / (sum_score_3 + sum_score_not_3),2)

    return mh_score, oof_score



'''
Test
target = 'GTGCGGCTGGCCCAGGACCTAGG'
target30mer = 'CTTGTGCGGCTGGCCCAGGACCTAGGCGAG'
target60mer = 'ATCTCTTACAACGACTTCTTGTGCGGCTGGCCCAGGACCTAGGCGAGGCAGTAGGGGATGACA'

print('Microhomology, Out-Of-Frame score: ',oofscore(seq = target60mer))
## Ans: 3730.1, 64.2

print('Azimuth on-target score: ',azimuth([target30mer])[0])
#Ans: 0.38
'''

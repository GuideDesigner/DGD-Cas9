import RNA
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict, OrderedDict
from Bio.SeqUtils import MeltingTemp as mt
import itertools
from itertools import chain
import make_arrays
import stacking_model
import string
import tensorflow as tf
from tensorflow.python.client import device_lib
import tensorflow.keras.backend as kb
from tensorflow.keras import models, layers, optimizers, losses

def test_model(modelD, X, a):
    resultList = []
    for modelF in filter(lambda x: x.endswith('h5'), os.listdir(modelD)):
        model = models.load_model(modelD + '/' + modelF)

        test_result = model.predict([X, a])

        resultList.append(test_result)

    # resultList = np.asarray(resultList)
    result = np.mean(resultList, axis=0)

    return (result)


def get_Input(inputF):
    trans = str.maketrans('ACGT', '0123')
    # writer = open(inputF.split('.csv')[0] + '_Sequence.txt','w')
    in_ = []
    add_ = []
    idss = []
    totalN = 0
    f = open(inputF)
    q = open(inputF)
    klines = q.readlines()
    lines = f.readlines()
    del klines[0]
    f.close()
    for line in klines:
        info_id = line.strip().split(',')[0]
        idss.append(info_id)

    colnames = [x.strip('"') for x in lines[0].strip().split(',')]
    Index_As = [colnames.index('A'+str(x)) for x in range(1, 31)]
    Index_Ts = [colnames.index('T'+str(x)) for x in range(1, 31)]
    Index_Gs = [colnames.index('G'+str(x)) for x in range(1, 31)]
    Index_Cs = [colnames.index('C'+str(x)) for x in range(1, 31)]

    if 'AA1' in colnames:
        Index_AAs = [colnames.index('AA'+str(x)) for x in range(1, 30)]
        Index_ATs = [colnames.index('AT'+str(x)) for x in range(1, 30)]
        Index_AGs = [colnames.index('AG'+str(x)) for x in range(1, 30)]
        Index_ACs = [colnames.index('AC'+str(x)) for x in range(1, 30)]
        Index_TAs = [colnames.index('TA'+str(x)) for x in range(1, 30)]
        Index_TTs = [colnames.index('TT'+str(x)) for x in range(1, 30)]
        Index_TGs = [colnames.index('TG'+str(x)) for x in range(1, 30)]
        Index_TCs = [colnames.index('TC'+str(x)) for x in range(1, 30)]
        Index_GAs = [colnames.index('GA'+str(x)) for x in range(1, 30)]
        Index_GTs = [colnames.index('GT'+str(x)) for x in range(1, 30)]
        Index_GGs = [colnames.index('GG'+str(x)) for x in range(1, 30)]
        Index_GCs = [colnames.index('GC'+str(x)) for x in range(1, 30)]
        Index_CAs = [colnames.index('CA'+str(x)) for x in range(1, 30)]
        Index_CTs = [colnames.index('CT'+str(x)) for x in range(1, 30)]
        Index_CGs = [colnames.index('CG'+str(x)) for x in range(1, 30)]
        Index_CCs = [colnames.index('CC'+str(x)) for x in range(1, 30)]
    else:
        Index_AAs = []
        Index_ATs = []
        Index_AGs = []
        Index_ACs = []
        Index_TAs = []
        Index_TTs = []
        Index_TGs = []
        Index_TCs = []
        Index_GAs = []
        Index_GTs = []
        Index_GGs = []
        Index_GCs = []
        Index_CAs = []
        Index_CTs = []
        Index_CGs = []
        Index_CCs = []

    Index_ID = [colnames.index('ID')]
    Index_others = [x for x in range(len(colnames)) if not x in Index_As + Index_Ts + Index_Gs + Index_Cs + Index_ID + Index_AAs+Index_ATs +
                    Index_AGs+Index_ACs+Index_TAs+Index_TTs+Index_TGs+Index_TCs+Index_GAs+Index_GTs+Index_GGs+Index_GCs+Index_CAs+Index_CTs+Index_CGs+Index_CCs]
    
    for line in lines[1:]:
        line = line.strip().split(',')
        _A = [int(line[x]) for x in Index_As]
        _T = [int(line[x]) for x in Index_Ts]
        _G = [int(line[x]) for x in Index_Gs]
        _C = [int(line[x]) for x in Index_Cs]
        _others = [float(line[x]) for x in Index_others]
        seq = ''
        for i in range(len(_A)):
            if _A[i] == 1:
                seq += 'A'
            elif _T[i] == 1:
                seq += 'T'
            elif _C[i] == 1:
                seq += 'C'
            elif _G[i] == 1:
                seq += 'G'
            else:
                print("Not in ATGC !")
        var_onehot = []
        for i in range(len(_A)):
            var_onehot.append([_A[i], _C[i], _G[i], _T[i]])
        kmerList = [seq[i:i+2] for i in range(len(seq)-1)]

        dinucList = []
        for nuc in itertools.product('ATCG', repeat=2):
            dinuc = ''.join(nuc)
            for kmer in kmerList:
                if kmer == dinuc:
                    dinucList.append(1)
                else:
                    dinucList.append(0)
        _others += dinucList

        in_.append(var_onehot)
        add_.append(_others)
    in_ = np.asarray(in_)
    in_ = in_.reshape(-1, 30, 4, 1)
    add_ = np.asarray(add_)
    ids = np.asarray(Index_ID)

    return in_, add_, idss


def score_deep():
    struct_file = pd.read_csv("Structure_file.csv")
    in_, add_, pop_id = get_Input(
        "Deep_learning_file.csv")

    SpCas9 = test_model('./SpCas9/', in_, add_)
    SpCas9_score = pd.DataFrame(SpCas9.reshape(len(SpCas9), 1))
    SpCas9_score.columns = ['DGDSpCas9']
    ID = pd.DataFrame(pop_id, columns=['ID'])
    Score_complete = pd.concat([ID, SpCas9_score], axis=1, join='outer')
    Merge_str = pd.merge(struct_file, Score_complete, on='ID', how='inner')
    Merge_str.to_csv("DGD.csv", index=False)

import os
import sys
import random
import itertools
import numpy as np
import tensorflow as tf
from tensorflow.python.client import device_lib
import tensorflow.keras.backend as kb
from tensorflow.keras import models, layers, optimizers, losses
import pandas as pd

def test_model(model, X, a):

    if len(a) > 0 and len(X) > 0:
        test_result = model.predict([X, a])
    elif len(X) == 0:
        test_result = model.predict(a)
    else:
        test_result = model.predict(X)
    return test_result


def score_deep(file1,file2):
    struct_file = pd.read_csv(file1)
    trans = str.maketrans('ACGT', '0123')
    #writer = open(inputF.split('.csv')[0] + '_Sequence.txt','w')
    idq = []
    in_ = []
    out_ = []
    add_ = []
    class_ = []
    totalN = 0
    f = open(file2, 'r')
    lines = f.readlines()
    f.close()
    colnames = [x.strip('"') for x in lines[0].strip().split(',')]
    Index_As = [colnames.index('A'+str(x)) for x in range(1, 31)]
    Index_Ts = [colnames.index('T'+str(x)) for x in range(1, 31)]
    Index_Gs = [colnames.index('G'+str(x)) for x in range(1, 31)]
    Index_Cs = [colnames.index('C'+str(x)) for x in range(1, 31)]
    Index_ID = [colnames.index('ID')]

    try:
        Index_indel = [colnames.index('INDEL')]
    except ValueError:
        Index_indel = None
    try:
        Index_class = [colnames.index('Class')]
    except ValueError:
        Index_class = None
    try:
        CC_class = [colnames.index('Connection_tab')]
    except ValueError:
        CC_class = None
    if Index_class == None and CC_class == None and Index_indel == None:
        Index_others = [x for x in range(
            len(colnames)) if not x in Index_As + Index_Ts + Index_Gs + Index_Cs + Index_ID]
    elif Index_class == None:
        Index_others = [x for x in range(len(colnames)) if not x in Index_As +
                        Index_Ts + Index_Gs + Index_Cs + Index_ID + Index_indel + CC_class]
    elif CC_class == None:
        Index_others = [x for x in range(len(colnames)) if not x in Index_As +
                        Index_Ts + Index_Gs + Index_Cs + Index_ID + Index_indel + Index_class]
    else:
        Index_others = [x for x in range(len(colnames)) if not x in Index_As + Index_Ts +
                        Index_Gs + Index_Cs + Index_ID + Index_indel + Index_class + CC_class]
    for line in lines[1:]:
        totalN += 1
        if totalN == 1:
            continue
        line = line.strip().split(',')
        ids = [(line[x]) for x in Index_ID][0]

        idq.append(ids)
 #       score = [float(line[x]) for x in Index_indel][0]
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

        class_.append(None)
        in_.append(var_onehot)
        out_.append([None])
        if len(_others) == 0:
            pass
        else:
            add_.append(_others)
    in_ = np.asarray(in_)
    in_ = in_.reshape(-1, 30, 4, 1)
#    out_ = np.asarray(out_).reshape(-1,1)
#    class_ = np.asarray(class_).reshape(-1,1)
    add_ = np.asarray(add_)
    idss = np.asarray(idq)
    
    global add_shape
    add_shape = add_.shape[1:]
    m1 = models.load_model('./models/Best_1fold.h5')
    m2 = models.load_model('./models/Best_2fold.h5')
    m3 = models.load_model('./models/Best_3fold.h5')
    m4 = models.load_model('./models/Best_4fold.h5')
    m5 = models.load_model('./models/Best_5fold.h5')
    m6 = models.load_model('./models/Best_6fold.h5')
    m7 = models.load_model('./models/Best_7fold.h5')
    m8 = models.load_model('./models/Best_8fold.h5')
    m9 = models.load_model('./models/Best_9fold.h5')
    m10 = models.load_model('./models/Best_10fold.h5')
    s1 = test_model(m1, in_, add_)
    s2 = test_model(m2, in_, add_)
    s3 = test_model(m3, in_, add_)
    s4 = test_model(m4, in_, add_)
    s5 = test_model(m5, in_, add_)
    s6 = test_model(m6, in_, add_)
    s7 = test_model(m7, in_, add_)
    s8 = test_model(m8, in_, add_)
    s9 = test_model(m9, in_, add_)
    s10 = test_model(m10, in_, add_)
    score_data_1 = pd.DataFrame(s1.reshape(len(s1), 1))
    score_data_2 = pd.DataFrame(s1.reshape(len(s2), 1))
    score_data_3 = pd.DataFrame(s1.reshape(len(s3), 1))
    score_data_4 = pd.DataFrame(s1.reshape(len(s4), 1))
    score_data_5 = pd.DataFrame(s1.reshape(len(s5), 1))
    score_data_6 = pd.DataFrame(s1.reshape(len(s6), 1))
    score_data_7 = pd.DataFrame(s1.reshape(len(s7), 1))
    score_data_8 = pd.DataFrame(s1.reshape(len(s8), 1))
    score_data_9 = pd.DataFrame(s1.reshape(len(s9), 1))
    score_data_10 = pd.DataFrame(s1.reshape(len(s10), 1))
    
    dgd_id = pd.DataFrame(idss.reshape(len(idss), 1))
    score_stack = pd.DataFrame(np.column_stack([score_data_1, score_data_2, score_data_3, score_data_4, score_data_5, score_data_6, score_data_7, score_data_8, score_data_9, score_data_10]), columns=[
                                   "Model_1", "Model_2", "Model_3", "Model_4", "Model_5", "Model_6", "Model_7", "Model_8", "Model_9", "Model_10"])
    score_stack['Score'] = score_stack.mean(axis=1)
    final_score_deep = pd.DataFrame(score_stack, columns=['Score'])
    final_score_stack = pd.DataFrame(np.column_stack([final_score_deep, dgd_id]), columns=["Score_Deep", "ID"])
    Merge_str = pd.merge(struct_file,final_score_stack,on='ID',how='inner')
    Merge_str.to_csv("DGD_score.csv", index=False)

args1 = sys.argv[1]
args2 = sys.argv[2]
result = score_deep(args1,args2)
print(result)

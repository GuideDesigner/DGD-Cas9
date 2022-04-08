################################################################################################# 
# Copyright (c),A Vipin Menon, Jang-il Sohn, Seokju Park and BIG LAB in Hanyang University(HYU) #
# Author A Vipin Menon 									        #
# Date 21st August, 2021								        #	
# Email a.vipin.menon@gmail.com				 				        #
#################################################################################################

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


def reverseString(st):
    li = []
    for i in st:
        li.append(i)
    li.reverse()
    return ''.join(li)


def reverseComp(st):
    comp = str.maketrans('ATCG', 'TAGC')
    return reverseString(st).translate(comp)


def targetsequence():
    sequence = open("Structure_file.csv", 'r')
    out = open("Target_sequence_feature.csv", 'w')
    header = list(map(lambda x: 'A' + str(x), range(1, 31))) + list(map(lambda x: 'T' + str(x), range(1, 31))) + list(map(lambda x: 'G' + str(x), range(1, 31))) + list(map(lambda x: 'C' + str(x), range(1, 31))) + list(map(lambda x: 'AA'+str(x), range(1, 30))) + list(map(lambda x: 'TA' + str(x), range(1, 30))) + list(map(lambda x: 'GA' + str(x), range(1, 30))) + list(map(lambda x: 'CA' + str(x), range(1, 30))) + list(map(lambda x: 'AT' + str(x), range(1, 30))) + list(map(lambda x: 'TT' + str(x), range(1, 30))) + \
        list(map(lambda x: 'GT' + str(x), range(1, 30))) + list(map(lambda x: 'CT' + str(x), range(1, 30))) + list(map(lambda x: 'AG' + str(x), range(1, 30))) + list(map(lambda x: 'TG' + str(x), range(1, 30))) + list(map(lambda x: 'GG' + str(x), range(1, 30))) + \
        list(map(lambda x: 'CG' + str(x), range(1, 30))) + list(map(lambda x: 'AC' + str(x), range(1, 30))) + list(map(lambda x: 'TC' +
                                                                                                                       str(x), range(1, 30))) + list(map(lambda x: 'GC' + str(x), range(1, 30))) + list(map(lambda x: 'CC' + str(x), range(1, 30)))
    out.write('ID' + ',')
    out.write(','.join(header))
    out.write(',' + 'Entropy' + ',' + 'Energy' + ',' + 'GCcount' + ',' + 'Gchigh' + ',' + 'GClow' + ',' + 'MeltingTemperature' + ',' + 'A' + ',' + 'T' + ',' + 'G' + ',' + 'C' + ',' + 'AA' + ',' + 'AT' +
              ',' + 'AG' + ',' + 'AC' + ',' + 'CA' + ',' + 'CG' + ',' + 'CC' + ',' + 'CT' + ',' + 'GA' + ',' + 'GC' + ',' + 'GG' + ',' + 'GT' + ',' + 'TA' + ',' + 'TC' + ',' + 'TG' + ',' + 'TT' + '\n')
    sequence = sequence.readlines()
    del sequence[0]
    l = []
    dinuct = []
    ene = {}
    ent = {}
    ext = {}
    dt = {}
    nt = {}
    PosA = {}
    PosT = {}
    PosC = {}
    PosG = {}
    Tm = {}
    PosAA = {}
    PosAT = {}
    PosAG = {}
    PosAC = {}
    PosCA = {}
    PosCC = {}
    PosCG = {}
    PosCT = {}
    PosGA = {}
    PosGC = {}
    PosGG = {}
    PosGT = {}
    PosTA = {}
    PosTC = {}
    PosTG = {}
    PosTT = {}
    gc = {}
    gc_high = {}
    gc_low = {}
    nuc = ['A', 'T', 'G', 'C']
    nas = []
    merged_dict = defaultdict(list)
    sigma_dict = defaultdict(list)
    # ntscount = {'A':0, 'G':0, 'T':0, 'C':0}
    for j in range(0, len(nuc)):
        for t in range(0, len(nuc)):
            p = str(nuc[t]) + str(nuc[j])
            dinuct.append(p)
    for line in sequence:
        seq = line.strip().split(',')
    #	number = seq[2]
        indel = seq[2]
        complete_sequence = seq[4]
        ids = seq[0]
        target_sequence = complete_sequence[4:24].upper()
    #	start = complete_sequence.find(target_sequence)
    #	newstart = start - 21
    #	newend = int(len(target_sequence))  + 21
        Extended_sequence = complete_sequence.upper()

    #	print len(Extended_sequence),Extended_sequence
        ext[ids] = indel
        nt[ids] = []
        dt[ids] = []
        ent[ids] = []
        ene[ids] = []
        PosA[ids] = []
        PosT[ids] = []
        PosC[ids] = []
        PosG[ids] = []
        PosAA[ids] = []
        PosAT[ids] = []
        PosAG[ids] = []
        PosAC[ids] = []
        PosGA[ids] = []
        PosGC[ids] = []
        PosGG[ids] = []
        PosGT[ids] = []
        PosCA[ids] = []
        PosCG[ids] = []
        PosCC[ids] = {}
        PosCT[ids] = []
        PosTA[ids] = []
        PosTG[ids] = []
        PosTC[ids] = []
        PosTT[ids] = []
        Tm[ids] = []
        gc[ids] = []
        gc_high[ids] = []
        gc_low[ids] = []
        for y in range(0, len(nuc)):
            for x in range(0, len(Extended_sequence)):
                if nuc[y] in Extended_sequence[x]:
                    Count = str(1)
                    nt[ids].append(Count)
                else:

                    Count = str(0)
                    nt[ids].append(Count)

        for w in range(0, len(dinuct)):
            for z in range(0, len(Extended_sequence)-1):
                if dinuct[w] == Extended_sequence[z:z+2]:
                    Count = str(1)
                    dt[ids].append(Count)
                else:
                    Count = str(0)
                    dt[ids].append(Count)

        entropy = dict()
        Entropycal = []
        entropy_seq = target_sequence
        lentseq = len(entropy_seq)
        ntscount = {'A': 0, 'G': 0, 'T': 0, 'C': 0}
        for ant in nuc:
            ntscount[ant] = (entropy_seq.count(ant))/float((lentseq))
        for ant in nuc:
            if ntscount[ant] != 0:
                entropy[ant] = -(ntscount[ant]*np.log2(ntscount[ant]))
            else:
                entropy[ant] = 0

        entropySum = sum(entropy.values())
        entropySumR = round(entropySum, 1)
        ent[ids] = str(entropySumR)
        Energy = target_sequence
        Energycal = RNA.fold(Energy)[-1]
        Energycal = round(Energycal, 0)
        ene[ids] = str(Energycal)
        PosA[ids] = str(Extended_sequence.count('A'))
        PosC[ids] = str(Extended_sequence.count('C'))
        PosT[ids] = str(Extended_sequence.count('T'))
        PosG[ids] = str(Extended_sequence.count('G'))
        Tm[ids] = str(mt.Tm_NN(target_sequence))
        PosA[ids] = str(Extended_sequence.count('A'))
        PosT[ids] = str(Extended_sequence.count('T'))
        PosG[ids] = str(Extended_sequence.count('G'))
        PosC[ids] = str(Extended_sequence.count('C'))
        PosAA[ids] = str(Extended_sequence.count('AA'))
        PosAT[ids] = str(Extended_sequence.count('AT'))
        PosAG[ids] = str(Extended_sequence.count('AG'))
        PosAC[ids] = str(Extended_sequence.count('AC'))
        PosCA[ids] = str(Extended_sequence.count('CA'))
        PosCC[ids] = str(Extended_sequence.count('CC'))
        PosCG[ids] = str(Extended_sequence.count('CG'))
        PosCT[ids] = str(Extended_sequence.count('CT'))
        PosCC[ids] = str(Extended_sequence.count('CC'))
        PosGA[ids] = str(Extended_sequence.count('GA'))
        PosGC[ids] = str(Extended_sequence.count('GC'))
        PosGG[ids] = str(Extended_sequence.count('GG'))
        PosGT[ids] = str(Extended_sequence.count('GT'))
        PosTA[ids] = str(Extended_sequence.count('TA'))
        PosTC[ids] = str(Extended_sequence.count('TC'))
        PosTG[ids] = str(Extended_sequence.count('TG'))
        PosTT[ids] = str(Extended_sequence.count('TT'))
        gc_content = (target_sequence.count(
            'G') + target_sequence.count('C'))/float(len(target_sequence)) * 100
        gc_content = round(gc_content, 0)
        gc_count = (target_sequence.count('G') + target_sequence.count('C'))
        gc[ids] = str(gc_content)
        if gc_count < int(10):
            alpha = str(1)
        else:
            alpha = str(0)
        if gc_count >= int(10):
            beta = str(1)
        else:
            beta = str(0)

        gc_low[ids] = alpha
        gc_high[ids] = beta

    # for t,v in nt.items():
    #	v = ''.join[v]
    for key, value in nt.items():
        nt[key] = ','.join(value)
    for key, value in dt.items():
        dt[key] = ','.join(value)
#	print nt,dt

    dict_list = [nt, dt, ent, ene, gc, gc_high, gc_low, Tm, PosA, PosT, PosG, PosC, PosAA, PosAT, PosAG,
                 PosAC, PosCA, PosCC, PosCG, PosCT, PosGA, PosGC, PosGG, PosGT, PosTA, PosTC, PosTG, PosTT]

    for dicts in dict_list:
        for k, v in dicts.items():
            merged_dict[k].append(v)
#	print merged_dict
    for key, value in (merged_dict.items()):
        sigma_dict[key] = ','.join((value))

    for key, value in sorted(sigma_dict.items()):
        out.write(str(key) + ',' + str(value) + '\n')
    out.close()


def Connectstr():
    f1 = open("Structure_Connection.outs", 'r')
    f1 = f1.readlines()
    mac = []
    sara = {}
    a = []
    kt = []
    header = map(lambda x: 'Pos' + str(x), range(1, 103))
    n = 'ID'
    a.append(n)
    a.extend(header)

    for line in f1:
        info = line.strip().split(' ')
        info_list = list(filter(None, info))

        if len(info_list) == 5:

            mykey = info_list[4]
            sara[mykey] = []

        else:
            sara[mykey].append(info_list[4])

    df = pd.DataFrame.from_dict(sara, orient='index')
    df.to_csv("Structure_Cas9_out.txt", sep="\t", header=False)
    df = pd.read_csv("Structure_Cas9_out.txt", sep="\t", names=a)
    df.to_csv("Structure_out.txt", sep="\t", index=False)
    os.remove("Structure_Cas9_out.txt")


def spacerscaffold():
    conn_d = pd.read_csv('Structure_basepairs.csv')
    a = []

    for i in range(1, 21):
        for j in range(21, 103):
            n = 'Connection_Pos' + str(i) + '_Pos' + str(j)
            a.append(n)
    a.extend(['ID'])

    newdf = pd.DataFrame(conn_d, columns=a)
    newdf.to_csv("spacer_scaffold_basepairs.csv", index=False)


def spacerconnectionfrequency():

    connection = pd.read_csv("spacer_scaffold_basepairs.csv")
    cons = connection.set_index('ID').T.reset_index()

    cons_info = cons[cons.columns[0]]
    cons_info = pd.DataFrame(cons_info)
    cons_info[['nucleotide', 'Pos_A', 'Pos_B']] = pd.DataFrame(
        [x.split('_') for x in cons_info[cons_info.columns[0]].tolist()])
    cons_info = cons_info.drop('nucleotide', axis=1)
    cons_info['Pos_A'] = cons_info['Pos_A'].str.extract('(\d+)', expand=False)
    cons_info['Pos_B'] = cons_info['Pos_B'].str.extract('(\d+)', expand=False)
    cons_info.columns = ['nucleotide', 'Pos_A', 'Pos_B']
    cons_info.to_csv("spacer_scaffold_feature.csv", index=False)


def serialconnection():
    spacer_scaffold = pd.read_csv('spacer_scaffold_feature.csv')
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 33) & (
        spacer_scaffold.Pos_B <= 36), 'Structure'] = 'TL'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 54) & (
        spacer_scaffold.Pos_B <= 58), 'Structure'] = 'SL1'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 73) & (
        spacer_scaffold.Pos_B <= 76), 'Structure'] = 'SL2'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 88) & (
        spacer_scaffold.Pos_B <= 90), 'Structure'] = 'SL3'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 21) & (
        spacer_scaffold.Pos_B <= 32), 'Structure'] = 'R'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 37) & (
        spacer_scaffold.Pos_B <= 49), 'Structure'] = 'AR'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 63) & (
        spacer_scaffold.Pos_B <= 67), 'Structure'] = 'LR'
    spacer_scaffold.Structure.fillna('NS', inplace=True)
    spacer_scaffold.to_csv("Structural_annotation.csv", index=False)


def featuremaker():

    Connection_bp = pd.read_csv("spacer_scaffold_basepairs.csv")
    Sequence = pd.read_csv("Target_sequence_feature.csv")
    Serial_Conn = pd.read_csv("Structural_annotation.csv")
    Serial_conn_info_NS = Serial_Conn[(Serial_Conn.Structure == 'NS')][[
        "nucleotide", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_AR = Serial_Conn[(Serial_Conn.Structure == 'AR')][[
        "nucleotide", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_R = Serial_Conn[(Serial_Conn.Structure == 'R')][[
        "nucleotide", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_LR = Serial_Conn[(Serial_Conn.Structure == 'LR')][[
        "nucleotide", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_SL1 = Serial_Conn[(Serial_Conn.Structure == 'SL1')][[
        "nucleotide", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_SL2 = Serial_Conn[(Serial_Conn.Structure == 'SL2')][[
        "nucleotide", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_SL3 = Serial_Conn[(Serial_Conn.Structure == 'SL3')][[
        "nucleotide", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_TL = Serial_Conn[(Serial_Conn.Structure == 'TL')][[
        "nucleotide", "Pos_A", "Pos_B", "Structure"]]

    position_list = [x for x in range(1, 21)]

    PosA = Serial_conn_info_AR[["Pos_A"]]

    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)

    nested_AR = {}
    C_AR = {}
    for x in range(0, len(newlist)):
        AR = {}
        name = 'AR' + str(newlist[x])
        my_list = Serial_conn_info_AR[(Serial_conn_info_AR.Pos_A == newlist[x])][[
            "nucleotide"]]

        mylist = my_list.values.tolist()
        palist = list(set([item for sublist in mylist for item in sublist]))
        Connection_information = Connection_bp[palist]
        value = Connection_bp[['ID']]
        Connection_information.insert(0, 'ID', value)

        my_id_dict = Connection_information.set_index('ID').T.to_dict('list')
        for key, value in my_id_dict.items():
            for x in range(0, len(value)):
                if sum(value) > 0:
                    AR[key] = 1
                else:
                    AR[key] = 0
        nested_AR[name] = (AR)

    AR = pd.DataFrame(nested_AR)
    AR['Total_AR'] = AR.sum(axis=1)
    AR.index.name = 'ID'
    AR.reset_index(inplace=True)
   # total_AR.index.name = 'ID'
   # total_AR.reset_index(inplace=True)

    PosA = Serial_conn_info_NS[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)

    nested_NS = {}
    C_NS = {}
    for x in range(0, len(newlist)):
        NS = {}
        name = 'NS' + str(newlist[x])
        my_list = Serial_conn_info_NS[(Serial_conn_info_NS.Pos_A == newlist[x])][[
            "nucleotide"]]
        mylist = my_list.values.tolist()
        palist = list(set([item for sublist in mylist for item in sublist]))
        Connection_information = Connection_bp[palist]
        value = Connection_bp[['ID']]
        Connection_information.insert(0, 'ID', value)

        my_id_dict = Connection_information.set_index('ID').T.to_dict('list')
        for key, value in my_id_dict.items():
            for x in range(0, len(value)):
                if sum(value) > 0:
                    NS[key] = 1
                    C_NS[key] = sum(value)
                else:
                    NS[key] = 0
        nested_NS[name] = (NS)
    NS = pd.DataFrame(nested_NS)
    NS['Total_NS'] = NS.sum(axis=1)
    NS.index.name = 'ID'
    NS.reset_index(inplace=True)

    PosA = Serial_conn_info_R[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)

    nested_R = {}
    C_R = {}
    for x in range(0, len(newlist)):
        R = {}
        name = 'R' + str(newlist[x])
        my_list = Serial_conn_info_R[(Serial_conn_info_R.Pos_A == newlist[x])][[
            "nucleotide"]]
        mylist = my_list.values.tolist()
        palist = list(set([item for sublist in mylist for item in sublist]))
        Connection_information = Connection_bp[palist]
        value = Connection_bp[['ID']]
        Connection_information.insert(0, 'ID', value)

        my_id_dict = Connection_information.set_index('ID').T.to_dict('list')
        for key, value in my_id_dict.items():
            for x in range(0, len(value)):
                if sum(value) > 0:
                    R[key] = 1
                    C_R[key] = sum(value)
                else:
                    R[key] = 0
        nested_R[name] = (R)
    R = pd.DataFrame(nested_R)
    R['Total_R'] = R.sum(axis=1)
    R.index.name = 'ID'
    R.reset_index(inplace=True)

    PosA = Serial_conn_info_SL1[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_SL1 = {}
    C_SL1 = {}
    for x in range(0, len(newlist)):
        SL1 = {}
        name = 'SL1_' + str(newlist[x])
        my_list = Serial_conn_info_SL1[(
            Serial_conn_info_SL1.Pos_A == newlist[x])][["nucleotide"]]
        mylist = my_list.values.tolist()
        palist = list(set([item for sublist in mylist for item in sublist]))
        Connection_information = Connection_bp[palist]
        value = Connection_bp[['ID']]
        Connection_information.insert(0, 'ID', value)

        my_id_dict = Connection_information.set_index('ID').T.to_dict('list')
        for key, value in my_id_dict.items():
            for x in range(0, len(value)):
                if sum(value) > 0:
                    SL1[key] = 1
                    C_SL1[key] = sum(value)
                else:
                    SL1[key] = 0
        nested_SL1[name] = SL1
    SL1 = pd.DataFrame(nested_SL1)
    SL1['Total_SL1'] = SL1.sum(axis=1)
    SL1.index.name = 'ID'
    SL1.reset_index(inplace=True)

    PosA = Serial_conn_info_LR[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_LR = {}
    C_LR = {}
    for x in range(0, len(newlist)):
        LR = {}
        name = 'LR' + str(newlist[x])
        my_list = Serial_conn_info_LR[(Serial_conn_info_LR.Pos_A == newlist[x])][[
            "nucleotide"]]
        mylist = my_list.values.tolist()
        palist = list(set([item for sublist in mylist for item in sublist]))
        Connection_information = Connection_bp[palist]
        value = Connection_bp[['ID']]
        Connection_information.insert(0, 'ID', value)

        my_id_dict = Connection_information.set_index('ID').T.to_dict('list')
        for key, value in my_id_dict.items():
            for x in range(0, len(value)):
                if sum(value) > 0:
                    LR[key] = 1
                    C_LR[key] = sum(value)
                else:
                    LR[key] = 0
        nested_LR[name] = LR
    LR = pd.DataFrame(nested_LR)
    LR['Total_LR'] = LR.sum(axis=1)
    LR.index.name = 'ID'
    LR.reset_index(inplace=True)

    PosA = Serial_conn_info_SL2[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_SL2 = {}
    C_SL2 = {}
    for x in range(0, len(newlist)):
        SL2 = {}
        name = 'SL2_' + str(newlist[x])
        my_list = Serial_conn_info_SL2[(
            Serial_conn_info_SL2.Pos_A == newlist[x])][["nucleotide"]]
        mylist = my_list.values.tolist()
        palist = list(set([item for sublist in mylist for item in sublist]))
        Connection_information = Connection_bp[palist]
        value = Connection_bp[['ID']]
        Connection_information.insert(0, 'ID', value)

        my_id_dict = Connection_information.set_index('ID').T.to_dict('list')
        for key, value in my_id_dict.items():
            for x in range(0, len(value)):
                if sum(value) > 0:
                    SL2[key] = 1
                    C_SL2[key] = sum(value)
                else:
                    SL2[key] = 0
        nested_SL2[name] = SL2
    SL2 = pd.DataFrame(nested_SL2)
    SL2['Total_SL2'] = SL2.sum(axis=1)
    SL2.index.name = 'ID'
    SL2.reset_index(inplace=True)

    PosA = Serial_conn_info_SL3[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_SL3 = {}
    C_SL3 = {}
    for x in range(0, len(newlist)):
        SL3 = {}
        name = 'SL3_' + str(newlist[x])
        my_list = Serial_conn_info_SL3[(
            Serial_conn_info_SL3.Pos_A == newlist[x])][["nucleotide"]]
        mylist = my_list.values.tolist()
        palist = list(set([item for sublist in mylist for item in sublist]))
        Connection_information = Connection_bp[palist]
        value = Connection_bp[['ID']]
        Connection_information.insert(0, 'ID', value)

        my_id_dict = Connection_information.set_index('ID').T.to_dict('list')
        for key, value in my_id_dict.items():
            for x in range(0, len(value)):
                if sum(value) > 0:
                    SL3[key] = 1
                    C_SL3[key] = sum(value)
                else:
                    SL3[key] = 0
        nested_SL3[name] = SL3
    SL3 = pd.DataFrame(nested_SL3)
    SL3['Total_SL3'] = SL3.sum(axis=1)
    SL3.index.name = 'ID'
    SL3.reset_index(inplace=True)

    PosA = Serial_conn_info_TL[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_TL = {}
    C_TL = {}
    for x in range(0, len(newlist)):
        TL = {}
        name = 'TL_' + str(newlist[x])
        my_list = Serial_conn_info_TL[(Serial_conn_info_TL.Pos_A == newlist[x])][[
            "nucleotide"]]
        mylist = my_list.values.tolist()
        palist = list(set([item for sublist in mylist for item in sublist]))
        Connection_information = Connection_bp[palist]
        value = Connection_bp[['ID']]
        Connection_information.insert(0, 'ID', value)

        my_id_dict = Connection_information.set_index('ID').T.to_dict('list')
        for key, value in my_id_dict.items():
            for x in range(0, len(value)):
                if sum(value) > 0:
                    TL[key] = 1
                    C_TL[key] = sum(value)
                else:
                    TL[key] = 0
        nested_TL[name] = TL
    TL = pd.DataFrame(nested_TL)
    TL['Total_TL'] = TL.sum(axis=1)
    TL.index.name = 'ID'
    TL.reset_index(inplace=True)
    Conn = pd.concat([R['Total_R'], TL['Total_TL'], AR['Total_AR'], SL1['Total_SL1'],
                      LR['Total_LR'], SL2['Total_SL2'], SL3['Total_SL3'], NS['Total_NS']], axis=1, join='outer')
    total = pd.DataFrame(Conn.sum(axis=1))

    total.columns = ["Total_Connection"]

    R = R.drop(["ID"], axis=1)
    LR = LR.drop(["ID"], axis=1)
    SL1 = SL1.drop(["ID"], axis=1)
    SL2 = SL2.drop(["ID"], axis=1)
    SL3 = SL3.drop(["ID"], axis=1)
    TL = TL.drop(["ID"], axis=1)
    NS = NS.drop(["ID"], axis=1)
    Connection = pd.concat(
        [AR, R, LR, SL2, TL, NS, SL1, SL3, total], axis=1, join='outer')
#    Connection[['ID', 'INDEL']] = Connection['ID'].str.split(':', 1, expand=True)
#    Connection['INDEL'] = Connection['INDEL'].astype(float)
#    Sequence['INDEL'] = Sequence['INDEL'].astype(float)
    Connection_sequence = Sequence.merge(Connection, on=['ID'], how='inner')

    Connection_sequence.to_csv("Feature_Data_Spacer_Scaffold.csv", index=False)


def read_dot_bracket_format_file(file_name):
    """
    Read a file in the RNAfold' Dot-bracket format:
    read 3 rows at each time, the first row is the sequence ID, the second row is the sequence, the third row is the Dot-bracket structure
    {'ID': [sequence, structure, free-energey], ...}
    """
    file = open(file_name, 'r')
    ID_sequence_structure_energy = {}
    lines = file.readlines()
    for i in range(len(lines)):
        if (i % 3 == 0):
            tmp = lines[i]
#            ID = tmp[1:tmp.find(':')]
            ID = tmp[1:].strip('\n')
            ID_sequence_structure_energy[ID] = []
        elif (i % 3 == 1):
            sequence = lines[i]
            ID_sequence_structure_energy[ID].append(
                sequence[:sequence.find('\n')])
        elif (i % 3 == 2):
            structure = lines[i]
            ID_sequence_structure_energy[ID].append(
                structure[:structure.rfind('(')])
            ID_sequence_structure_energy[ID].append(
                structure[structure.rfind('('):structure.find('\n')])
    file.close()
    return ID_sequence_structure_energy


# count mono-nucleotide
def count_monomer(sequence, partner_index):
    """
    Input: sequence, partner_base
    Return the number of mono-nucleotide
    """
    monomer = {}
    monomer = {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    for i in range(20):
        if partner_index[i] != 0:
            monomer[sequence[i]] += 1
            monomer[sequence[partner_index[i]]] += 1
    return monomer


# count di-nucleotide
def count_dimer(dimers):
    """
    Input: sequence, partner_base
    Return the number of di-nucleotide
    """
    dinucleotide = {}
    dinucleotide = {'AA': 0, 'AC': 0, 'AG': 0, 'AU': 0,
                    'CA': 0, 'CC': 0, 'CG': 0, 'CU': 0,
                    'GA': 0, 'GC': 0, 'GG': 0, 'GU': 0,
                    'UA': 0, 'UC': 0, 'UG': 0, 'UU': 0}

    for i in range(len(dimers)):
        dinucleotide[dimers[i]] += 1
    return dinucleotide


# count GC content
def count_GC_content(sequence, partner_index):
    """
    Count GC content of a sequence
    """
    count = 0
    for i in range(20):
        if partner_index[i] != 0:
            if (sequence[i] == 'G' and sequence[partner_index[i]] == 'C'):
                count += 1
            elif (sequence[i] == 'C' and sequence[partner_index[i]] == 'G'):
                count += 1
    return count


def calculate_free_energy(dimers):
    """
    Calculate free energy of a dimer, return the free energy
    First, convert the sequence and partner_base to dimers arrays
    Second, calculate the free energy of each dimer based on stacking model
    """
    free_energy = 0

    stacking_model_array = stacking_model.return_stacking_model()
    for i in range(0, len(dimers), 2):
        free_energy += stacking_model_array[dimers[i]][dimers[i+1][::-1]]
    return free_energy


def finalfeatures():
    ID_sequence_structure_energy = read_dot_bracket_format_file(
        "Structure_Connection.out")
    feature_data = pd.read_csv("Feature_Data_Spacer_Scaffold.csv")

    IDS = {}
    IDX = {}
    monomer = []
    dimer = []
    for ID in ID_sequence_structure_energy:
        sequence = ID_sequence_structure_energy[ID][0]
        structure = ID_sequence_structure_energy[ID][1]
        energy = ID_sequence_structure_energy[ID][2]

        partner_index = make_arrays.return_partner_index(structure)  # num
        partner_index = make_arrays.adjust_partner_index(partner_index)
        partner_base = make_arrays.return_partner_base(
            partner_index, sequence)  # seq
        connection_length, count = make_arrays.return_connection_length(
            partner_index)

        # make dimers array
        dimers = make_arrays.return_dimers_array(sequence, partner_index)

        # mask array
        SP_sequence_masked = make_arrays.return_masked_array(
            sequence, partner_index, 20)  # num
        SP_partner_base_masked = make_arrays.return_masked_array(
            partner_base, partner_index, 20)  # num

        # count GC content
        if count != 0:
            GC_count = count_GC_content(sequence, partner_index)/count
        else:
            GC_count = 0
        # count mono-nucleotide
        monomer_count = count_monomer(sequence, partner_index)
        for key, values in monomer_count.items():
            monomer.append(values)

        # count dimers
        dimers_count = count_dimer(dimers)
        for key, values in dimers_count.items():
            dimer.append(values)

        # stacking energy
        stacking_energy = calculate_free_energy(dimers)
        IDS[ID] = [GC_count, stacking_energy]
        IDX[ID] = [monomer_count, dimers_count]

        monomer_dimer = {k: {k: v for IDX in L for k, v in IDX.items()}
                         for k, L in IDX.items()}

    df_1 = pd.DataFrame(monomer_dimer).T
    df_1.index.name = 'ID'
    df_1.reset_index(inplace=True)

    df_2 = pd.DataFrame.from_dict(IDS).T
    df_2.columns = ["GC_ratio", "Gibbs_Energy"]
    df_2.index.name = 'ID'
    df_2.reset_index(inplace=True)

    mydf = pd.merge(df_1, df_2, on="ID", how='inner')
    mydf.columns = ['ID', 'Spacer_Scaffold_A', 'Spacer_Scaffold_C', 'Spacer_Scaffold_G', 'Spacer_Scaffold_U', 'Spacer_Scaffold_AA', 'Spacer_Scaffold_AC', 'Spacer_Scaffold_AG', 'Spacer_Scaffold_AU', 'Spacer_Scaffold_CA', 'Spacer_Scaffold_CC', 'Spacer_Scaffold_CG',
                    'Spacer_Scaffold_CU', 'Spacer_Scaffold_GA', 'Spacer_Scaffold_GC', 'Spacer_Scaffold_GG', 'Spacer_Scaffold_GU', 'Spacer_Scaffold_UA', 'Spacer_Scaffold_UC', 'Spacer_Scaffold_UG', 'Spacer_Scaffold_UU', 'Spacer_Scaffold_GC_ratio', 'Spacer_Scaffold_Gibbs_Energy']
#    feature_data['INDEL'] = feature_data['INDEL'].astype(float)

    final_feature = feature_data.merge(mydf, on=['ID'])
    final_feature.to_csv("Deep_learning_file.csv", index=False)


def fastamaker():
    file1 = open("Structure_file.csv", 'r')
    out = open('Structure_Connection.fa', 'w')
    gout = open('Structure_Connection.csv', 'w')
    gout.write('ID' + ',' + 'Sequence' + '\n')
    file1 = file1.readlines()
    seq = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT'
    del file1[0]
    for line in file1:
        info = line.strip().split(',')
        ids = info[0]
        sequence = info[4][4:24] + seq
        out.write('>' + str(ids) + '\n' + str(sequence) + '\n')
        gout.write(str(ids) + ',' + str(sequence) + '\n')
    out.close()
    gout.close()


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


def dgdmain(file1):
    cas9 = {}
    ncas9 = {}
    sara = {}
    ast = []
    fastafile = open(file1, 'r')
    input1 = file1[:-3]
    out = open("Structure_file.csv", 'w')
    fastafile = fastafile.readlines()
    out.write('ID' + ',' + 'Start' + ',' + 'End' +
              ',' + 'Strand' + ',' + 'Sequence' + '\n')
    for line in fastafile:
        if line.startswith(">"):
            newline = line.strip('\r\n')
            nare = newline[1:]
            nare = nare.split(' ')[0]
            ps = nare
            sara[ps] = []
        else:
            newsequence = line.strip('\r\n').split('\n')
            for az in newsequence:
                ast.append(az)
            my_str = ''.join(map(str, ast))
            sara[ps] = [my_str]
    for key, value in sara.items():
        my_value = value[0]
        if (int(100) <= len(my_value) <= int(10000)):
            for n in range(len(my_value)):
                st = my_value.find('GG', n)
                if st == n:
                    if int(25) < int(st) < len(my_value)-int(5):
                        nstart = int(st) - int(25)
                        nend = int(st) + int(5)
                        sequences = my_value[nstart:nend]
                        idn = key + ':' + str(nstart) + ':' + str(nend)
                        cas9[idn] = [nstart, nend, '+', sequences]
            for p in range(len(my_value)):
                st = my_value.find('CC', p)
                if st == p:
                    if int(3) < int(st) < len(my_value)-int(27):
                        nstart = int(st) - int(3)
                        nend = int(st) + int(27)
                        sequences = my_value[nstart:nend]
                        nsequences = reverseComp(sequences)
                        idn = key + ':' + str(nstart) + ':' + str(nend)
                        ncas9[idn] = [nstart, nend, '-', nsequences]
    newcasi = dict(chain(cas9.items(), ncas9.items()))

    syt = OrderedDict(
        sorted(newcasi.items(), key=lambda x: x[1][3], reverse=True))

    for key, value in syt.items():
        out.write(str(
            key) + ',' + str(value[0]) + ',' + str(value[1]) + ',' + str(value[2]) + ',' + str(value[3]) + '\n')
    out.close()

    fastamaker()
    targetsequence()
    os.system(
        "RNAfold -j0 --noPS <Structure_Connection.fa> Structure_Connection.out")
    os.system("b2ct <Structure_Connection.out> Structure_Connection.outs")
    Connectstr()
    os.system(
        "./connection_to_matrix Structure_out.txt 102 > Structure_basepairs.csv")
    spacerscaffold()
    spacerconnectionfrequency()
    serialconnection()
    featuremaker()
    finalfeatures()
    score_deep()


args = sys.argv[1]
result = dgdmain(args)
print(result)

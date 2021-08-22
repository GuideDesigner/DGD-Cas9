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
from itertools import chain
import string
import tensorflow as tf
from tensorflow.python.client import device_lib
import tensorflow.keras.backend as kb
from tensorflow.keras import models, layers, optimizers, losses



def reverseString(st):
        li = []
        for i in st: li.append(i)
        li.reverse()
        return ''.join(li)


def reverseComp(st):
        comp = str.maketrans('ATCG', 'TAGC')
        return reverseString(st).translate(comp)


def fastamaker():
	file1 = open("Structure_file.csv",'r')
	out = open('Structure_Connection.fa','w')
	gout = open('Structure_Connection.csv','w')
	gout.write('ID' + ',' + 'Sequence' + '\n')
	file1 = file1.readlines()
	seq = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT'
	del file1[0]
	for line in file1:
		info = line.strip().split(',')
		ids = info[0]
		sequence = info[4][4:24]  + seq
		out.write('>' + str(ids) + '\n' + str(sequence) + '\n')
		gout.write(str(ids) + ',' + str(sequence) + '\n')
	out.close()
	gout.close()
def sequence():
    sequence = open("Structure_file.csv", 'r')
    out = open("Target_sequence_feature.csv", 'w')
    header = list(map(lambda x: 'A' + str(x), range(1, 31))) + list(map(lambda x: 'T' + str(x), range(1, 31))) + list(map(lambda x: 'G' + str(x), range(1, 31))) + list(map(lambda x: 'C' + str(x), range(1, 31))) + list(map(lambda x: 'AA'+str(x), range(1, 30))) + list(map(lambda x: 'TA'+ str(x), range(1, 30))) + list(map(lambda x: 'GA' + str(x), range(1, 30))) + list(map(lambda x: 'CA' + str(x), range(1, 30))) + list(map(lambda x: 'AT' + str(x), range(1, 30))) + list(map(lambda x: 'TT' + str(x), range(1, 30))) + list(map(lambda x: 'GT' + str(x), range(1, 30))) + list(map(lambda x: 'CT' + str(x), range(1, 30))) + list(map(lambda x: 'AG' + str(x), range(1, 30))) + list(map(lambda x: 'TG' + str(x), range(1, 30))) + list(map(lambda x: 'GG' + str(x), range(1, 30))) + list(map(lambda x: 'CG' + str(x), range(1, 30))) + list(map(lambda x: 'AC' + str(x), range(1, 30))) + list(map(lambda x: 'TC' + str(x), range(1, 30))) + list(map(lambda x: 'GC' + str(x), range(1, 30))) + list(map(lambda x: 'CC' + str(x), range(1, 30)))
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
#        indel = seq[2]
        complete_sequence = seq[4]
        ids = seq[0]
        target_sequence = complete_sequence[4:24]
    #	start = complete_sequence.find(target_sequence)
    #	newstart = start - 21
    #	newend = int(len(target_sequence))  + 21
        Extended_sequence = complete_sequence
    #	print len(Extended_sequence),Extended_sequence
    #    ext[ids] = indel
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
        Tm[ids] = str(float(64.9) + float(41.0) * ((float(target_sequence.count('G')
                                                          ) + float(target_sequence.count('C')) - float(16.4))/(20.0)))
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
    df = pd.read_csv("Structure_Cas9_out.txt", sep="\t")
    df.columns = a
    df.to_csv("Structure_out.txt", sep="\t", index=False)
    os.remove("Structure_Cas9_out.txt")

def spacerscaffold ():
	conn_d = pd.read_csv('Structure_basepairs.csv')	
	a = []

	for i in range(1,21):
		for j in range(21,103):
			n = 'Connection_Pos' + str(i) + '_Pos' + str(j)
			a.append(n)
	a.extend(['ID'])

	newdf = pd.DataFrame(conn_d,columns = a)
	newdf.to_csv("spacer_scaffold_basepairs.csv",index=False)

def spacerconnectionfrequency():

	connection = pd.read_csv("spacer_scaffold_basepairs.csv")
	cons = connection.set_index('ID').T.reset_index()

	cons_info = cons[cons.columns[0]]
	cons_info = pd.DataFrame(cons_info)


	cons_info[['nucleotide','Pos_A','Pos_B']] = pd.DataFrame([x.split('_') for x in cons_info[cons_info.columns[0]].tolist()])
	cons_info = cons_info.drop('nucleotide',axis=1)
	cons_info['Pos_A'] = cons_info['Pos_A'].str.extract('(\d+)',expand=False)
	cons_info['Pos_B'] = cons_info['Pos_B'].str.extract('(\d+)',expand=False)
	cons_info.columns = ['nucleotide', 'Pos_A', 'Pos_B']
	cons_info.to_csv("spacer_scaffold_feature.csv",index=False)

def featuremaker():

    Connection_bp = pd.read_csv("spacer_scaffold_basepairs.csv")
    Sequence = pd.read_csv("Target_sequence_feature.csv")
    Serial_Conn = pd.read_csv("Serial_connection_spacer_scaffold.csv")
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
    AR.index.name = 'ID'
    AR.reset_index(inplace=True)

    PosA = Serial_conn_info_NS[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)

    nested_NS = {}
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
                else:
                    NS[key] = 0
        nested_NS[name] = (NS)
    NS = pd.DataFrame(nested_NS)
    NS.index.name = 'ID'
    NS.reset_index(inplace=True)

    PosA = Serial_conn_info_R[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)

    nested_R = {}
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
                else:
                    R[key] = 0
        nested_R[name] = (R)
    R = pd.DataFrame(nested_R)
    R.index.name = 'ID'
    R.reset_index(inplace=True)

    PosA = Serial_conn_info_SL1[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_SL1 = {}
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
                else:
                    SL1[key] = 0
        nested_SL1[name] = SL1
    SL1 = pd.DataFrame(nested_SL1)
    SL1.index.name = 'ID'
    SL1.reset_index(inplace=True)

    PosA = Serial_conn_info_LR[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_LR = {}
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
                else:
                    LR[key] = 0
        nested_LR[name] = LR
    LR = pd.DataFrame(nested_LR)
    LR.index.name = 'ID'
    LR.reset_index(inplace=True)

    PosA = Serial_conn_info_SL2[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_SL2 = {}
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
                else:
                    SL2[key] = 0
        nested_SL2[name] = SL2
    SL2 = pd.DataFrame(nested_SL2)
    SL2.index.name = 'ID'
    SL2.reset_index(inplace=True)

    PosA = Serial_conn_info_SL3[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_SL3 = {}
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
                else:
                    SL3[key] = 0
        nested_SL3[name] = SL3
    SL3 = pd.DataFrame(nested_SL3)
    SL3.index.name = 'ID'
    SL3.reset_index(inplace=True)

    PosA = Serial_conn_info_TL[["Pos_A"]]
    PosA = PosA.values.tolist()
    flat_list = list(set([item for sublist in PosA for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    nested_TL = {}
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
                else:
                    TL[key] = 0
        nested_TL[name] = TL
    TL = pd.DataFrame(nested_TL)
    TL.index.name = 'ID'
    TL.reset_index(inplace=True)

    unique_id = Serial_conn_info_AR[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    CC = {}
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[['ID']]

    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_AR = pd.DataFrame.from_dict(CC, orient='index')
    name = 'AR'
    CC_AR.columns = [name]
    CC_AR.index.name = 'ID'
    CC_AR.reset_index(inplace=True)

    unique_id = Serial_conn_info_R[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))
    newlist = sorted(flat_list, reverse=True)
    CC = {}
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[['ID']]

    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_R = pd.DataFrame.from_dict(CC, orient='index')
    name = 'R'
    CC_R.columns = [name]
    CC_R.index.name = 'ID'
    CC_R.reset_index(inplace=True)

    unique_id = Serial_conn_info_LR[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))
    newlist = sorted(flat_list, reverse=True)
    CC = {}
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[['ID']]
    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_LR = pd.DataFrame.from_dict(CC, orient='index')
    name = 'LR'
    CC_LR.columns = [name]
    CC_LR.index.name = 'ID'
    CC_LR.reset_index(inplace=True)

    unique_id = Serial_conn_info_TL[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))

    newlist = sorted(flat_list, reverse=True)
    CC = {}
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[['ID']]
    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_TL = pd.DataFrame.from_dict(CC, orient='index')
    name = 'TL'
    CC_TL.columns = [name]
    CC_TL.index.name = 'ID'
    CC_TL.reset_index(inplace=True)

    unique_id = Serial_conn_info_SL2[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))
    newlist = sorted(flat_list, reverse=True)
    CC = {}
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[['ID']]
    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_SL2 = pd.DataFrame.from_dict(CC, orient='index')
    name = 'SL2'
    CC_SL2.columns = [name]
    CC_SL2.index.name = 'ID'
    CC_SL2.reset_index(inplace=True)

    unique_id = Serial_conn_info_SL1[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))
    newlist = sorted(flat_list, reverse=True)
    CC = {}
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[['ID']]
    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_SL1 = pd.DataFrame.from_dict(CC, orient='index')
    name = 'SL1'
    CC_SL1.columns = [name]
    CC_SL1.index.name = 'ID'
    CC_SL1.reset_index(inplace=True)

    unique_id = Serial_conn_info_SL3[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))
    newlist = sorted(flat_list, reverse=True)
    CC = {}
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[['ID']]
    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_SL3 = pd.DataFrame.from_dict(CC, orient='index')
    name = 'SL3'
    CC_SL3.columns = [name]
    CC_SL3.index.name = 'ID'
    CC_SL3.reset_index(inplace=True)

    unique_id = Serial_conn_info_NS[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))
    newlist = sorted(flat_list, reverse=True)
    CC = {}
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[['ID']]
    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_NS = pd.DataFrame.from_dict(CC, orient='index')
    name = 'NS'
    CC_NS.columns = [name]
    CC_NS.index.name = 'ID'
    CC_NS.reset_index(inplace=True)

    unique_id = Serial_Conn[["nucleotide"]]
    unique_id = unique_id.values.tolist()
    flat_list = list(set([item for sublist in unique_id for item in sublist]))
    newlist = sorted(flat_list, reverse=True)
    CC = {}
    name = 'Connection_all'
    Connection_bp_information = Connection_bp[newlist]
    value = Connection_bp[["ID"]]
    Connection_bp_information.insert(0, 'ID', value)
    my_id_dict = Connection_bp_information.set_index('ID').T.to_dict('list')
    for key, value in my_id_dict.items():
        for x in range(0, len(value)):
            CC[key] = sum(value)

    CC_all = pd.DataFrame.from_dict(CC, orient='index')
    CC_all.columns = [name]
    CC_all.index.name = 'ID'
    CC_all.reset_index(inplace=True)

    CC_LR = CC_LR.drop(["ID"], axis=1)
    CC_R = CC_R.drop(["ID"], axis=1)
    CC_TL = CC_TL.drop(["ID"], axis=1)
    CC_NS = CC_NS.drop(["ID"], axis=1)
    CC_SL2 = CC_SL2.drop(["ID"], axis=1)
    CC_SL1 = CC_SL1.drop(["ID"], axis=1)
    CC_SL3 = CC_SL3.drop(["ID"], axis=1)
    CC_all = CC_all.drop(["ID"], axis=1)

    AR = AR.drop(["ID"], axis=1)
    R = R.drop(["ID"], axis=1)
    LR = LR.drop(["ID"], axis=1)
    SL1 = SL1.drop(["ID"], axis=1)
    SL2 = SL2.drop(["ID"], axis=1)
    SL3 = SL3.drop(["ID"], axis=1)
    TL = TL.drop(["ID"], axis=1)
    NS = NS.drop(["ID"], axis=1)
    CC_id_merged = pd.concat([CC_AR, CC_R, CC_LR, CC_SL2, CC_TL, CC_NS, CC_SL1,
                              CC_SL3, AR, R, LR, SL1, SL2, SL3, TL, NS, CC_all], axis=1, join='outer')

    CC_id_deep = pd.merge(Sequence, CC_id_merged, on='ID', how='inner')

    CC_id_deep.to_csv("Feature_Data_Spacer_Scaffold.csv", index=False)

def Finalfeatures():


	Sequence_tracr = pd.read_csv("Structure_Connection.csv")
	Serial_Conn = pd.read_csv("Serial_connection_spacer_scaffold.csv")
	Deep_learning = pd.read_csv("Feature_Data_Spacer_Scaffold.csv")

	Serial_Conn_Cons = Serial_Conn[(Serial_Conn.CC_num != 1)][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_TL = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'TL')][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_SL1 = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'SL1')][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_SL2 = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'SL2')][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_SL3 = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'SL3')][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_R = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'R')][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_AR = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'AR')][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_LR = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'LR')][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_NS = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'NS')][["nucleotide", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	energy_5 = {'AA': '-4.26', 'AT': '-3.67', 'TA': '-2.50', 'CA': '-6.12', 'GT': '-6.09','CT': '-5.40', 'GA': '-5.51', 'CG': '-9.07', 'GC': '-9.36', 'GG': '-7.66'}
	energy_3 = {'TT': '-4.26', 'TA': '-3.67', 'AT': '-2.50', 'GT': '-6.12', 'CA': '-6.09','GA': '-5.40', 'CT': '-5.51', 'GC': '-9.07', 'CG': '-9.36', 'CC': '-7.66'}
	term_5 = {'A': '4.3', 'G': '4.05'}
	term_3 = {'T': '4.3', 'C': '4.05'}

	unique_id = Serial_Conn_Cons_TL[["Unique_ID"]]
	unique_id = unique_id.values.tolist()
	flat_list = list(set([item for sublist in unique_id for item in sublist]))
	spacr_tracr = Sequence_tracr.set_index('ID').T.to_dict('list')

	newlist = sorted(flat_list, reverse=True)
	nestedspac = {}
	nestedA = {}
	nestedT = {}
	nestedG = {}
	nestedC = {}
	nestedAA = {}
	nestedAT = {}
	nestedAG = {}
	nestedAC = {}
	nestedGA = {}
	nestedGG = {}
	nestedGC = {}
	nestedGT = {}
	nestedTA = {}
	nestedTG = {}
	nestedTC = {}
	nestedTT = {}
	nestedCA = {}
	nestedCC = {}
	nestedCG = {}
	nestedCT = {}
	nestedGCcount = {}

	for x in range(0, len(newlist)):
		spac = {}
		A = {}
		T = {}
		G = {}
		C = {}
		AA = {}
		AG = {}
		AT = {}
		AC = {}
		GA = {}
		GG = {}
		GT = {}
		GC = {}
		CA = {}
		CG = {}
		CT = {}
		CC = {}
		TA = {}
		TG = {}
		TT = {}
		TC = {}
		GCC = {}
		A_name = "A_" + str(newlist[x])
		T_name = "T_" + str(newlist[x])
		G_name = "G_" + str(newlist[x])
		C_name = "C_" + str(newlist[x])
		AA_name = "AA_" + str(newlist[x])
		AG_name = "AG_" + str(newlist[x])
		AT_name = "AT_" + str(newlist[x])
		AC_name = "AC_" + str(newlist[x])
		TA_name = "TA_" + str(newlist[x])
		TG_name = "TG_" + str(newlist[x])
		TC_name = "TC_" + str(newlist[x])
		TT_name = "TT_" + str(newlist[x])
		GA_name = "GA_" + str(newlist[x])
		GG_name = "GG_" + str(newlist[x])
		GT_name = "GT_" + str(newlist[x])
		GC_name = "GC_" + str(newlist[x])
		CA_name = "CA_" + str(newlist[x])
		CG_name = "CG_" + str(newlist[x])
		CT_name = "CT_" + str(newlist[x])
		CC_name = "CC_" + str(newlist[x])
		GC_count = "GC_count_" + str(newlist[x])
		name = 'CC_' + str(newlist[x])

		PosA = Serial_Conn_Cons[(
		Serial_Conn_Cons.Unique_ID == newlist[x])][["Pos_A"]]
		posa = PosA.values.tolist()
		PosB = Serial_Conn_Cons[(
		Serial_Conn_Cons.Unique_ID == newlist[x])][["Pos_B"]]
		posb = PosB.values.tolist()
		posa_list = list(set([item for sublist in posa for item in sublist]))
		posb_list = list(set([item for sublist in posb for item in sublist]))
		spacr_start = int(min(posa_list)) - int(1)
		spacr_end = int(max(posa_list))
		tracr_start = int(min(posb_list)) - int(1)
		tracr_end = int(max(posb_list))

		for key, value in spacr_tracr.items():
			spacrseq = value[0][spacr_start:spacr_end]
			tracrseq = value[0][tracr_start:tracr_end]
			sac = str(spacrseq) + str(tracrseq)
			gc = int(sac.count('G')) + int(sac.count('C'))
			rnafold = RNA.fold(sac)[-1]
			spac[key] = rnafold
			A[key] = (sac.count('A'))
			T[key] = (sac.count('T'))
			G[key] = (sac.count('G'))
			C[key] = (sac.count('C'))
			AA[key] = (sac.count('AA'))
			AT[key] = (sac.count('AT'))
			AG[key] = (sac.count('AG'))
			AC[key] = (sac.count('AC'))
			GA[key] = (sac.count('GA'))
			GT[key] = (sac.count('GT'))
			GC[key] = (sac.count('GC'))
			GG[key] = (sac.count('GG'))
			TA[key] = (sac.count('TA'))
			TT[key] = (sac.count('TT'))
			TC[key] = (sac.count('TC'))
			TG[key] = (sac.count('TG'))
			CA[key] = (sac.count('CA'))
			CT[key] = (sac.count('CT'))
			CC[key] = (sac.count('CC'))
			CG[key] = (sac.count('CG'))
			GCC[key] = gc
		nestedspac[name] = spac
		nestedA[A_name] = A
		nestedT[T_name] = T
		nestedG[G_name] = G
		nestedC[C_name] = C
		nestedAA[AA_name] = AA
		nestedAT[AT_name] = AT
		nestedAG[AG_name] = AG
		nestedAC[AC_name] = AC
		nestedGA[GA_name] = GA
		nestedGG[GG_name] = GG
		nestedGC[GC_name] = GC
		nestedGT[GT_name] = GT
		nestedTA[TA_name] = TA
		nestedTC[TC_name] = TC
		nestedTT[TT_name] = TT
		nestedTG[TG_name] = TG
		nestedCA[CA_name] = CA
		nestedCC[CC_name] = CC
		nestedCG[CG_name] = CG
		nestedCT[CT_name] = CT
		nestedGCcount[GC_count] = GCC

	Energy = pd.DataFrame(nestedspac)
	Energy.index.name = 'ID'
	Energy.reset_index(inplace=True)

	Energy_ID = Energy[["ID"]]
	Energy_rest = Energy.drop(columns=['ID'])
	Energy_rest = Energy_rest.mean(axis=1)
	Energy_rest = pd.DataFrame(Energy_rest)
	Energy_rest.columns = ["Spacer_Scaffold_MFE"]
	value = Energy_ID[["ID"]]
	Energy_rest.insert(0, 'ID', value)

	NucA = pd.DataFrame(nestedA)
	NucT = pd.DataFrame(nestedT)
	NucG = pd.DataFrame(nestedG)
	NucC = pd.DataFrame(nestedC)
	NucAA = pd.DataFrame(nestedAA)
	NucAT = pd.DataFrame(nestedAT)
	NucAG = pd.DataFrame(nestedAG)
	NucAC = pd.DataFrame(nestedAC)
	NucTA = pd.DataFrame(nestedTA)
	NucTT = pd.DataFrame(nestedTT)
	NucTG = pd.DataFrame(nestedTG)
	NucTC = pd.DataFrame(nestedTC)
	NucCA = pd.DataFrame(nestedCA)
	NucCT = pd.DataFrame(nestedCT)
	NucCG = pd.DataFrame(nestedCG)
	NucCC = pd.DataFrame(nestedCC)
	NucGA = pd.DataFrame(nestedGA)
	NucGT = pd.DataFrame(nestedGT)
	NucGG = pd.DataFrame(nestedGG)
	NucGC = pd.DataFrame(nestedGC)
	gccount = pd.DataFrame(nestedGCcount)

	NucA.index.name = 'ID'
	NucA.reset_index(inplace=True)
	NucA_ID = NucA[["ID"]]
	NucA_rest = NucA.drop(columns=["ID"])
	NucA_rest = NucA_rest.mean(axis=1)
	NucA_rest = pd.DataFrame(NucA_rest)
	NucA_rest.columns = ["Spacer_Scaffold_A"]

	NucT.index.name = 'ID'
	NucT.reset_index(inplace=True)
	NucT_ID = NucT[["ID"]]
	NucT_rest = NucT.drop(columns=["ID"])
	NucT_rest = NucT_rest.mean(axis=1)
	NucT_rest = pd.DataFrame(NucT_rest)
	NucT_rest.columns = ["Total_T"]

	NucG.index.name = 'ID'
	NucG.reset_index(inplace=True)
	NucG_ID = NucG[["ID"]]
	NucG_rest = NucG.drop(columns=["ID"])
	NucG_rest = NucG_rest.mean(axis=1)
	NucG_rest = pd.DataFrame(NucG_rest)
	NucG_rest.columns = ["Spacer_Scaffold_G"]

	NucC.index.name = 'ID'
	NucC.reset_index(inplace=True)
	NucC_ID = NucC[["ID"]]
	NucC_rest = NucC.drop(columns=["ID"])
	NucC_rest = NucC_rest.mean(axis=1)
	NucC_rest = pd.DataFrame(NucC_rest)
	NucC_rest.columns = ["Spacer_Scaffold_C"]

	NucAA.index.name = 'ID'
	NucAA.reset_index(inplace=True)
	NucAA_ID = NucAA[["ID"]]
	NucAA_rest = NucAA.drop(columns=["ID"])
	NucAA_rest = NucAA_rest.mean(axis=1)
	NucAA_rest = pd.DataFrame(NucAA_rest)
	NucAA_rest.columns = ["Spacer_Scaffold_AA"]

	NucAT.index.name = 'ID'
	NucAT.reset_index(inplace=True)
	NucAT_ID = NucAT[["ID"]]
	NucAT_rest = NucAT.drop(columns=["ID"])
	NucAT_rest = NucAT_rest.mean(axis=1)
	NucAT_rest = pd.DataFrame(NucAT_rest)
	NucAT_rest.columns = ["Spacer_Scaffold_AT"]

	NucAG.index.name = 'ID'
	NucAG.reset_index(inplace=True)
	NucAG_ID = NucAG[["ID"]]
	NucAG_rest = NucAG.drop(columns=["ID"])
	NucAG_rest = NucAG_rest.mean(axis=1)
	NucAG_rest = pd.DataFrame(NucAG_rest)
	NucAG_rest.columns = ["Spacer_Scaffold_AG"]

	NucAC.index.name = 'ID'
	NucAC.reset_index(inplace=True)
	NucAC_ID = NucAC[["ID"]]
	NucAC_rest = NucAC.drop(columns=["ID"])
	NucAC_rest = NucAC_rest.mean(axis=1)
	NucAC_rest = pd.DataFrame(NucAC_rest)
	NucAC_rest.columns = ["Spacer_Scaffold_AC"]

	NucCA.index.name = 'ID'
	NucCA.reset_index(inplace=True)
	NucCA_ID = NucCA[["ID"]]
	NucCA_rest = NucCA.drop(columns=["ID"])
	NucCA_rest = NucCA_rest.mean(axis=1)
	NucCA_rest = pd.DataFrame(NucCA_rest)
	NucCA_rest.columns = ["Spacer_Scaffold_CA"]

	NucCG.index.name = 'ID'
	NucCG.reset_index(inplace=True)
	NucCG_ID = NucCG[["ID"]]
	NucCG_rest = NucCG.drop(columns=["ID"])
	NucCG_rest = NucCG_rest.mean(axis=1)
	NucCG_rest = pd.DataFrame(NucCG_rest)
	NucCG_rest.columns = ["Spacer_Scaffold_CG"]

	NucCT.index.name = 'ID'
	NucCT.reset_index(inplace=True)
	NucCT_ID = NucCT[["ID"]]
	NucCT_rest = NucCT.drop(columns=["ID"])
	NucCT_rest = NucCT_rest.mean(axis=1)
	NucCT_rest = pd.DataFrame(NucCT_rest)
	NucCT_rest.columns = ["Spacer_Scaffold_CT"]

	NucCC.index.name = 'ID'
	NucCC.reset_index(inplace=True)
	NucCC_ID = NucCC[["ID"]]
	NucCC_rest = NucCC.drop(columns=["ID"])
	NucCC_rest = NucCC_rest.mean(axis=1)
	NucCC_rest = pd.DataFrame(NucCC_rest)
	NucCC_rest.columns = ["Spacer_Scaffold_CC"]

	NucGG.index.name = 'ID'
	NucGG.reset_index(inplace=True)
	NucGG_ID = NucGG[["ID"]]
	NucGG_rest = NucGG.drop(columns=["ID"])
	NucGG_rest = NucGG_rest.mean(axis=1)
	NucGG_rest = pd.DataFrame(NucGG_rest)
	NucGG_rest.columns = ["Spacer_Scaffold_GG"]

	NucGT.index.name = 'ID'
	NucGT.reset_index(inplace=True)
	NucGT_ID = NucGT[["ID"]]
	NucGT_rest = NucGT.drop(columns=["ID"])
	NucGT_rest = NucGT_rest.mean(axis=1)
	NucGT_rest = pd.DataFrame(NucGT_rest)
	NucGT_rest.columns = ["Spacer_Scaffold_ _GT"]

	NucGA.index.name = 'ID'
	NucGA.reset_index(inplace=True)
	NucGA_ID = NucGA[["ID"]]
	NucGA_rest = NucGA.drop(columns=["ID"])
	NucGA_rest = NucGA_rest.mean(axis=1)
	NucGA_rest = pd.DataFrame(NucGA_rest)
	NucGA_rest.columns = ["Spacer_Scaffold_GA"]

	NucGC.index.name = 'ID'
	NucGC.reset_index(inplace=True)
	NucGC_ID = NucGC[["ID"]]
	NucGC_rest = NucGC.drop(columns=["ID"])
	NucGC_rest = NucGC_rest.mean(axis=1)
	NucGC_rest = pd.DataFrame(NucGC_rest)
	NucGC_rest.columns = ["Spacer_Scaffold_GC"]

	NucTC.index.name = 'ID'
	NucTC.reset_index(inplace=True)
	NucTC_ID = NucTC[["ID"]]
	NucTC_rest = NucTC.drop(columns=["ID"])
	NucTC_rest = NucTC_rest.mean(axis=1)
	NucTC_rest = pd.DataFrame(NucTC_rest)
	NucTC_rest.columns = ["Spacer_Scaffold_TC"]

	NucTA.index.name = 'ID'
	NucTA.reset_index(inplace=True)
	NucTA_ID = NucTA[["ID"]]
	NucTA_rest = NucTA.drop(columns=["ID"])
	NucTA_rest = NucTA_rest.mean(axis=1)
	NucTA_rest = pd.DataFrame(NucTA_rest)
	NucTA_rest.columns = ["Spacer_Scaffold_TA"]

	NucTG.index.name = 'ID'
	NucTG.reset_index(inplace=True)
	NucTG_ID = NucTG[["ID"]]
	NucTG_rest = NucTG.drop(columns=["ID"])
	NucTG_rest = NucTG_rest.mean(axis=1)
	NucTG_rest = pd.DataFrame(NucTG_rest)
	NucTG_rest.columns = ["Spacer_Scaffold_TG"]

	NucTT.index.name = 'ID'
	NucTT.reset_index(inplace=True)
	NucTT_ID = NucTT[["ID"]]
	NucTT_rest = NucTT.drop(columns=["ID"])
	NucTT_rest = NucTT_rest.mean(axis=1)
	NucTT_rest = pd.DataFrame(NucTT_rest)
	NucTT_rest.columns = ["Spacer_Scaffold_TT"]

	gccount.index.name = 'ID'
	gccount.reset_index(inplace=True)
	gccount_ID = gccount[["ID"]]
	gccount_rest = gccount.drop(columns=["ID"])
	gccount_rest = gccount_rest.mean(axis=1)
	gccount_rest = pd.DataFrame(gccount_rest)
	gccount_rest.columns = ["Spacer_Scaffold_GC_count"]

	CC_Total = pd.concat([Energy_rest, gccount_rest, NucA_rest, NucT_rest, NucG_rest, NucC_rest, NucAA_rest, NucAT_rest, NucAG_rest, NucAC_rest, NucGA_rest,
	  NucGT_rest, NucGG_rest, NucGC_rest, NucCA_rest, NucCT_rest, NucCG_rest, NucCC_rest, NucTA_rest, NucTT_rest, NucTG_rest, NucTC_rest], axis=1, join='outer')

	#print (CC_merse)

	CC_merge = pd.merge(CC_Total, Deep_learning, on='ID', how='inner')

	CC_merge.to_csv(
	"Deep_learning_feature_with_spacer_scaffold.csv", index=False)


def test_model(model, X, a):

    if len(a) > 0 and len(X) > 0:
        test_result = model.predict([X, a])
    elif len(X) == 0:
        test_result = model.predict(a)
    else:
        test_result = model.predict(X)
    return test_result


def score_deep():
    struct_file = pd.read_csv("Structure_file.csv")
    trans = str.maketrans('ACGT', '0123')
    #writer = open(inputF.split('.csv')[0] + '_Sequence.txt','w')
    idq = []
    in_ = []
    out_ = []
    add_ = []
    class_ = []
    totalN = 0
    f = open("Deep_learning_feature_with_spacer_scaffold.csv", 'r')
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
    final_score_stack = pd.DataFrame(np.column_stack([final_score_deep, dgd_id]), columns=["DGD", "ID"])
    Merge_str = pd.merge(struct_file,final_score_stack,on='ID',how='inner')
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
    out.write('ID' + ',' + 'Start' + ',' + 'End' + ',' + 'Strand' + ',' + 'Sequence' + '\n')
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
                        cas9[idn] = [nstart, nend, '+',sequences]
            for p in range(len(my_value)):
                st = my_value.find('CC', p)
                if st == p:
                    if int(3) < int(st) < len(my_value)-int(27):
                        nstart = int(st) - int(3)
                        nend = int(st) + int(27)
                        sequences = my_value[nstart:nend]
                        nsequences = reverseComp(sequences)
                        idn = key + ':' + str(nstart) + ':' + str(nend)
                        ncas9[idn] = [nstart, nend, '-',nsequences]
    newcasi = dict(chain(cas9.items(), ncas9.items()))

    syt = OrderedDict(
        sorted(newcasi.items(), key=lambda x: x[1][3], reverse=True))

    for key, value in syt.items():
        out.write(str(
            key) + ',' + str(value[0]) + ',' + str(value[1]) + ',' + str(value[2]) +  ',' + str(value[3]) + '\n')
    out.close()

    fastamaker()
    sequence()
    os.system("RNAfold -j0 --noPS <Structure_Connection.fa> Structure_Connection.out")
    os.system("b2ct <Structure_Connection.out> Structure_Connection.outs")
    Connectstr()   
    os.system("./connection_to_matrix Structure_out.txt 102 > Structure_basepairs.csv") 
    spacerscaffold()
    spacerconnectionfrequency()
    os.system("./CC_csv spacer_scaffold_feature.csv > CC_feature.csv")
    os.system("Rscript Serial_connection_training.R")
    featuremaker() 
    Finalfeatures() 
    score_deep()
    
args = sys.argv[1]
result = dgdmain(args)
print (result)

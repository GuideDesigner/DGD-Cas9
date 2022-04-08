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


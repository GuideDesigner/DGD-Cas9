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

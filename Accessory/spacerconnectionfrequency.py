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

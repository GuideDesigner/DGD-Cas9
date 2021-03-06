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

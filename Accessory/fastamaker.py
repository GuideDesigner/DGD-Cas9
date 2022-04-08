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

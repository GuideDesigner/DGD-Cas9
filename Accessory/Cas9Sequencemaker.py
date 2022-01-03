import sys
import RNA
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict, OrderedDict
from get_sequence import reverseString, reverseComp
from itertools import chain

def Cas9sequence(file1):
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

args = sys.argv[1]
result = Cas9sequence(args)

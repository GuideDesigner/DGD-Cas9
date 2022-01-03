import sys
import pandas as pd
import os


def Connectstr(file1):
    f1 = open(file1, 'r')
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

args = sys.argv[1]
result  =Connectstr(args)
print (result)


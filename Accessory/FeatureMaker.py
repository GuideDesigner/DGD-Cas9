import pandas as pd
import sys
import numpy as np


def featuremaker(x,y,z):

    Connection_bp = pd.read_csv("spacer_scaffold_basepairs.csv")
    Sequence = pd.read_csv("Target_sequence_feature.csv")
    Serial_Conn = pd.read_csv("Serial_connection_spacer_scaffold.csv")
    Serial_conn_info_NS = Serial_Conn[(Serial_Conn.Structure == 'NS')][[
        "connectionposition", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_AR = Serial_Conn[(Serial_Conn.Structure == 'AR')][[
        "connectionposition", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_R = Serial_Conn[(Serial_Conn.Structure == 'R')][[
        "connectionposition", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_LR = Serial_Conn[(Serial_Conn.Structure == 'LR')][[
        "connectionposition", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_SL1 = Serial_Conn[(Serial_Conn.Structure == 'SL1')][[
        "connectionposition", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_SL2 = Serial_Conn[(Serial_Conn.Structure == 'SL2')][[
        "connectionposition", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_SL3 = Serial_Conn[(Serial_Conn.Structure == 'SL3')][[
        "connectionposition", "Pos_A", "Pos_B", "Structure"]]
    Serial_conn_info_TL = Serial_Conn[(Serial_Conn.Structure == 'TL')][[
        "connectionposition", "Pos_A", "Pos_B", "Structure"]]

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
            "connectionposition"]]
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
            "connectionposition"]]
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
            "connectionposition"]]
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
            Serial_conn_info_SL1.Pos_A == newlist[x])][["connectionposition"]]
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
            "connectionposition"]]
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
            Serial_conn_info_SL2.Pos_A == newlist[x])][["connectionposition"]]
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
            Serial_conn_info_SL3.Pos_A == newlist[x])][["connectionposition"]]
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
            "connectionposition"]]
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

    unique_id = Serial_conn_info_AR[["connectionposition"]]
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

    unique_id = Serial_conn_info_R[["connectionposition"]]
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

    unique_id = Serial_conn_info_LR[["connectionposition"]]
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

    unique_id = Serial_conn_info_TL[["connectionposition"]]
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

    unique_id = Serial_conn_info_SL2[["connectionposition"]]
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

    unique_id = Serial_conn_info_SL1[["connectionposition"]]
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

    unique_id = Serial_conn_info_SL3[["connectionposition"]]
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

    unique_id = Serial_conn_info_NS[["connectionposition"]]
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

    unique_id = Serial_Conn[["connectionposition"]]
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

args1 = sys.argv[1]
args2 = sys.argv[2]
args3 = sys.argv[3]
result = featuremaker(args1,args2,args3)
print (result)

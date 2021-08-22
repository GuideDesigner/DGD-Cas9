import pandas as pd
import numpy as np
import RNA
import sys

def Finalfeatures(x,y,z):


	Sequence_tracr = pd.read_csv(x)
	Serial_Conn = pd.read_csv(y)
	Deep_learning = pd.read_csv(z)

	Serial_Conn_Cons = Serial_Conn[(Serial_Conn.CC_num != 1)][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_TL = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'TL')][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_SL1 = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'SL1')][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_SL2 = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'SL2')][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_SL3 = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'SL3')][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_R = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'R')][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_AR = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'AR')][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_LR = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'LR')][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

	Serial_Conn_Cons_NS = Serial_Conn_Cons[(Serial_Conn_Cons.Structure == 'NS')][["connectionposition", "Pos_A", "Pos_B", "CC_id", "CC_num", "Unique_ID", "Structure"]]

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

args1 = sys.argv[1]
args2 = sys.argv[2]
args3 = sys.argv[3]
result = Finalfeatures(args1,args2,args3)
print (result)

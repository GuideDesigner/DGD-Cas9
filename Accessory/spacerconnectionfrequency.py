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

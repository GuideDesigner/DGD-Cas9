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

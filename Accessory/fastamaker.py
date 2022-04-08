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

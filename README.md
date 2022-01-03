<<<<<<< HEAD
# DGD
=======
# DGD (Deep Guide Designer)
[![Build Status](https://travis-ci.org/joemccann/dillinger.svg?branch=master)](https://travis-ci.org/joemccann/dillinger)

DGD (Deep Guide Designer) is python based on-target scoring method for Cas9 sequences, which is developed by integrating information from target sequence, biofeatures and spacer-scaffold interaction of sgRNA into CNN based deep learning model.

## Source Code

### Briefings of python files
- Github directory contains main files and subsidary python files for the ease of users.
- Name of main file is ```DGD.py``` which when run will execute all the necessary files and give score as output.
- The input file should be fasta sequence (```.fa```), DGD pipeline can take multiple fasta sequence files.
- Accessory files (python) are provided for users in **Accessory** folder, who would like to study and integrate the features and outputs in their pipelines.
- Prior to use of any files both main and accessory, please install the packages as mentioned in ```Requirement``` section.

#### Requirements
- Python 3.6 or above
- Vienna RNA package
- Numpy (version >= 1.12.4)
- Pandas (version >= 1.1.4)
- Keras (version >= 2.2.5)
- Tensorflow(version == 2.2.0)
- R (version >= 3.6)
- gcc (version >= 4.7.0)

#### Package Installation
You should compile two c++ files:
```
make
```
Pre-compiled files(```connection_to_matrix``` and ```CC_csv```) are provided for ```DGD.py```

For Vienna RNA package see the instruction in their official site https://www.tbi.univie.ac.at/RNA/ for python installation
The rest packages could be installed by using ```pip``` or ```conda```

### Use
DGD python files is executed on fasta file
- The input file should in fasta format
- For demo purpose the length of fasta file is in range between 100nt to 10000nt
- The output are shown in ```.csv``` format which compromises of ID,Start,End,Strand, 30nt sequence and on-target score
```python
   python DGD.py input.fa
   
   Example file: input.fa
   
    >JQ277699.1 Ovis aries LIN28A (LIN28A) mRNA, complete cds
    ACCCCTTTGCCTCCGGACTTCTCTGGGGCCAGCAGCCGCCCAAGCGAGGGCCCGGGGCCGCGAGCTCAGC
    AGACTACCATGGGCTCTGTGTCAAACCAGCAGTTCGCAGGTGGCTGCGCTAAGGCGCCGGAGGAAGCGCC
    GGAGGAGGCGCCGGAGGACGCGGCCCGCGCGGCGGAGGAGCCGCAGCTGCTGCACGGTGCCGGCATCTGT
    AAGTGGTTCAACGTGCGCATGGGGTTCGGCTTCCTGTCCATGACCGCCCGCGCAGGGGTCGCGCTCGACC
    CCCCGGTGGATGTCTTTGTGCACCAGAGTAAGCTGCACATGGAGGGCTTCCGGAGCCTGAAGGAGGGGGA
    GGCTGTGGAGTTCACCTTTAAGAAGTCCGCCAAGGGCCTGGAATCTATCCGAGTCACCGGCCCTGGTGGG
    GTGTTCTGTATTGGGAGTGAAAGGCGGCCCAAAGGGAAGAATGTGCAGAAACGCAGATCAAAGGGAGACA
    GGTGCTACAACTGTGGAGGTCTAGACCATCATGCCAAGGAATGCAAACTGCCACCGCAGCCCAAGAAGTG
    CCATTTCTGCCAGAGCATCAGCCATATGGTAGCCTCGTGCCCACTGAAGGCCCAGCAAGCTCCCAGCTCC
    CAGGGAAAGCCAGCCTACTTTCGGGAGGAGGAAGAAGAGATCCATAGCTCTGCCATGCTCCCAGAGGCCC
    AGAATTGAAGCCACAGTGGGTGGGAGCTATCCTTTTGTGATCAGAAAGCTTTGAGGAGCAGCATCAATCG
    CCTTTGCCTTCGGACTTCTCCGGGGCCAGCAGCCGCCCGACCAGGGGCCCGGGGCCACGGGCTCAGCCGA
    CGACCATGGGCTCCGTGTCCAACCAGCAGTTTGCAGGTGGCTGCGCCAAGGCGGCAGAAGAGGCGCCCGA
    GGAGGCGCCGGAGGACGCGGCCCGGGCGGCGGACGAGCCTCAGCTGCTGCACGGTGCGGGCATCTGTAAG
    TGGTTCAACGTGCGCATGGGGTTCGGCTTCCTGTCCATGACCGCCCGCGCCGGGGTCGCGCTCGACCCCC
    CAGTGGATGTCTTTGTGCACCAGAGTAAGCTGCACATGGAAGGGTTCCGGAGCTTGAAGGAGGGTGAGGC
    AGTGGAGTTCACCTTTAAGAAGTCAGCCAAGGGTCTGGAATCCATCCGTGTCACCGGACCTGGTGGAGTA
    TTCTGTATTGGGAGTGAGAGGCGGCCAAAAGGAAAGAGCATGCAGAAGCGCAGATCAAAAGGAGACAGGT
    GCTACAACTGTGGAGGTCTAGATCATCATGCCAAGGAATGCAAGCTGCCACCCCAGCCCAAGAAGTGCCA

    >Q278699.2 Ovis aries 
    ACCCCTTTGCCTCCGGACTTCTCTGGGGCCAGCAGCCGCCCAAGCGAGGGCCCGGGGCCGCGAGCTCAGC
    AGACTACCATGGGCTCTGTGTCAAACCAGCAGTTCGCAGGTGGCTGCGCTAAGGCGCCGGAGGAAGCGCC
    GGAGGAGGCGCCGGAGGACGCGGCCCGCGCGGCGGAGGAGCCGCAGCTGCTGCACGGTGCCGGCATCTGT
    AAGTGGTTCAACGTGCGCATGGGGTTCGGCTTCCTGTCCATGACCGCCCGCGCAGGGGTCGCGCTCGACC
    CCCCGGTGGATGTCTTTGTGCACCAGAGTAAGCTGCACATGGAGGGCTTCCGGAGCCTGAAGGAGGGGGA
    GGCTGTGGAGTTCACCTTTAAGAAGTCCGCCAAGGGCCTGGAATCTATCCGAGTCACCGGCCCTGGTGGG
    GTGTTCTGTATTGGGAGTGAAAGGCGGCCCAAAGGGAAGAATGTGCAGAAACGCAGATCAAAGGGAGACA
    GGTGCTACAACTGTGGAGGTCTAGACCATCATGCCAAGGAATGCAAACTGCCACCGCAGCCCAAGAAGTG
    CCATTTCTGCCAGAGCATCAGCCATATGGTAGCCTCGTGCCCACTGAAGGCCCAGCAAGCTCCCAGCTCC
    CAGGGAAAGCCAGCCTACTTTCGGGAGGAGGAAGAAGAGATCCATAGCTCTGCCATGCTCCCAGAGGCCC
    AGAATTGAAGCCACAGTGGGTGGGAGCTATCCTTTTGTGATCAGAAAGCTTTGAGGAGCAGCATCAATCG
    CCTTTGCCTTCGGACTTCTCCGGGGCCAGCAGCCGCCCGACCAGGGGCCCGGGGCCACGGGCTCAGCCGA
```
- Output (DGD.csv)
```
ID	Start	End	Strand	Sequence	DGD
Q278699.2:1938:1968	1938	1968	-	TTTCCCTGGGAGCTGGGAGCTTGCTGGGCC	18.69614792
JQ277699.1:295:325	295	325	+	TTGTGCACCAGAGTAAGCTGCACATGGAGG	66.6506424
Q278699.2:295:325	295	325	+	TTGTGCACCAGAGTAAGCTGCACATGGAGG	66.6506424
Q278699.2:1625:1655	1625	1655	+	TTGTGCACCAGAGTAAGCTGCACATGGAGG	66.6506424
JQ277699.1:1062:1092	1062	1092	+	TTGTGCACCAGAGTAAGCTGCACATGGAAG	67.54428864
Q278699.2:1062:1092	1062	1092	+	TTGTGCACCAGAGTAAGCTGCACATGGAAG	67.54428864
JQ277699.1:859:889	859	889	-	TTGGCGCAGCCACCTGCAAACTGCTGGTTG	34.56448746
Q278699.2:859:889	859	889	-	TTGGCGCAGCCACCTGCAAACTGCTGGTTG	34.56448746
JQ277699.1:636:666	636	666	-	TTCTTCCTCCTCCCGAAAGTAGGCTGGCTT	47.40981674
Q278699.2:636:666	636	666	-	TTCTTCCTCCTCCCGAAAGTAGGCTGGCTT	47.40981674
Q278699.2:1966:1996	1966	1996	-	TTCTTCCTCCTCCCGAAAGTAGGCTGGCTT	47.40981674
JQ277699.1:1007:1037	1007	1037	+	TTCCTGTCCATGACCGCCCGCGCCGGGGTC	50.25638199
Q278699.2:1007:1037	1007	1037	+	TTCCTGTCCATGACCGCCCGCGCCGGGGTC	50.25638199
JQ277699.1:240:270	240	270	+	TTCCTGTCCATGACCGCCCGCGCAGGGGTC	60.52968216
Q278699.2:240:270	240	270	+	TTCCTGTCCATGACCGCCCGCGCAGGGGTC	60.52968216
Q278699.2:1570:1600	1570	1600	+	TTCCTGTCCATGACCGCCCGCGCAGGGGTC	60.52968216
JQ277699.1:607:637	607	637	-	TTCCCTGGGAGCTGGGAGCTTGCTGGGCCT	32.67510605
Q278699.2:607:637	607	637	-	TTCCCTGGGAGCTGGGAGCTTGCTGGGCCT	32.67510605
Q278699.2:1937:1967	1937	1967	-	TTCCCTGGGAGCTGGGAGCTTGCTGGGCCT	32.67510605
JQ277699.1:578:608	578	608	-	TTCAGTGGGCACGAGGCTACCATATGGCTG	59.14524078
Q278699.2:578:608	578	608	-	TTCAGTGGGCACGAGGCTACCATATGGCTG	59.14524078
Q278699.2:1908:1938	1908	1938	-	TTCAGTGGGCACGAGGCTACCATATGGCTG	59.14524078
JQ277699.1:679:709	679	709	-	TTCAATTCTGGGCCTCTGGGAGCATGGCAG	52.64225006
Q278699.2:679:709	679	709	-	TTCAATTCTGGGCCTCTGGGAGCATGGCAG	52.64225006
Q278699.2:2009:2039	2009	2039	-	TTCAATTCTGGGCCTCTGGGAGCATGGCAG	52.64225006
JQ277699.1:92:122	92	122	-	TTAGCGCAGCCACCTGCGAACTGCTGGTTT	19.55876541
Q278699.2:92:122	92	122	-	TTAGCGCAGCCACCTGCGAACTGCTGGTTT	19.55876541
Q278699.2:1422:1452	1422	1452	-	TTAGCGCAGCCACCTGCGAACTGCTGGTTT	19.55876541
JQ277699.1:86:116	86	116	+	TGTGTCAAACCAGCAGTTCGCAGGTGGCTG	55.41802216
```
## License
MIT
## Contributors
A. Vipin Menon, Jang-il Sohn, Seokju Park and Jin-Wu Nam

Bioinformatic and Genomics Lab., Hanyang University, Seoul 04763, Korea

## Contact
If you have any issues, please contact us

A. Vipin Menon (a.vipin.menon@gmail.com)

Jin-Wu Nam (jwnam@hanyang.ac.kr)
>>>>>>> 9f724f21e04c781da091c6c50276290ff64e40b9

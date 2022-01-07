'''
From a fasta file, find a specific transcript and its corresponding sequence.
usage: python <scriptName.py> <transcriptome.fasta> <ENSTid>
ex)TNFa, ENST00000374550.8
'''


import sys
infile=open(sys.argv[1])
id=open(sys.argv[2])


lines = infile.readlines()
for i in range(len(lines)):
    if lines[i].strip() == ">"+id:
        print (lines[i]+'\n'+lines[i + 1])

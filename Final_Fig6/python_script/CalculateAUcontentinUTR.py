'''AU content(%) in the 3' UTR list'''
'''usage: python <scriptName> <3UTR_concat.fasta> <outputName.txt>
outputFile format looks like as follows.
1st column : ENST Name
2nd column : A+U counts/ total length

'''

import sys
inFile=open(sys.argv[1])
outFile=open(sys.argv[2],"w")

for line in inFile:
    line=line.strip()
    if line[0]==">":
        id=line[1:]
        outFile.write(id+'\t')
    else:
        seq=line.strip()
        ATcontent=(seq.count('A')+seq.count('T'))/float(len(seq))
        outFile.write(str(ATcontent)+'\n')

outFile.close()

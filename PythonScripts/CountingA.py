'''counting A in a 6-mer list'''
'''
import sys

inFile=open(sys.argv[1])
outFile=open(sys.argv[2],"w")

for line in inFile:
    kmer=line.strip()
    ATcontent=(kmer.count('A')+kmer.count('T'))#/float(6)
    outFile.write(str(ATcontent)+'\n')
outFile.close()
'''



# AT contents in the 3'UTR list.

import sys
inFile=open(sys.argv[1])
outFile=open(sys.argv[2],"w")

for line in inFile:
    if line[0]==">":
        id=line[2:]
        outFile.write(id+'\t')
    else:
        seq=line.strip()
        ATcontent=seq.count('A')+seq.count('T')/float(len(seq))
        outFile.write(str(ATcontent)+'\n')

outFile.close()

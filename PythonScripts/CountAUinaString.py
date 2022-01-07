'''counting A in a 6-mer list'''
''' Usage: python <ScriptName.py> <6mer list.txt> <outputFileName>
ex) input: AATGCA ---> genearate output  3
Appliclation: This script can also be used to calculate AU content in any string file.
ex) fasta file:  every other line, AU count should be calculated and divided by length of string.


'''
import sys

inFile=open(sys.argv[1])
outFile=open(sys.argv[2],"w")

for line in inFile:
    kmer=line.strip()
    ATcontent=float(kmer.count('A')+kmer.count('T'))#/float(6)
    outFile.write(str(ATcontent)+'\n')
outFile.close()

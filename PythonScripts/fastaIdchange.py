'''Change too long fasta id format to concise one.'''
'''
Usage: python <scripName> <INPUTFILE> <OUTPUTFILE>
Action looks like as follows;
Change too long fasta id format to concise one.
input: >ENST00000332831.4_cds_0_0_chr1_685716_r::chr1:685715-686654(-)
output: >ENST00000332831.4
'''
import sys
inFile=open(sys.argv[1])
outFile=open(sys.argv[2],"w")

for line in inFile:
    line = line.strip()
    if line[0]== ">":
        l = line.split('_')
        outFile.write(l[0]+'\n')
    else:
        outFile.write(line+'\n')

outFile.close()

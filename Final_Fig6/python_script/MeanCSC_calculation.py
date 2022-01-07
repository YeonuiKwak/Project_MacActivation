'''usage:
python3 CAIanalysis_test.py ~/Desktop/fasta/GRCh38_CDS.fa ~/Desktop/ReadcountTable/Codon_CSC.txt CSC_scoreperid.txt'''


from itertools import chain
#import Bio.Data.CodonTable as ct
from scipy.stats import gmean
from collections import Counter
import sys

inFile=open(sys.argv[1])
#test=open(sys.argv[2])
#test=test.readlines()
sequences=inFile.readlines() # list
ref={}
for i in range(0,len(sequences),2):
    ref[sequences[i].strip()[1:]]=sequences[i+1].strip()

# count the number of each codon in the sequences
print(list(ref.items())[1:5])

no_three=0
ref2={}
#print ("total transcripts are:"+str(len(ref))) #106652
for key in ref:
    seq = ref[key].upper()
    if len(seq)%3==0 and ('N' not in seq):
        codon=[]
        for i in range(0,len(seq),3):
            codon.append(seq[i:i+3])
        ref2[key]=codon
    else:
        no_three+=1

print("total transcripts that have remainder when divided by 3 or has N in sequence:"+str(no_three)) #21422

#counts have fequencey of each codon in all transcipts.
#codons = chain.from_iterable(ref2) # flat list of all codons (to be used for counting)
#counts = Counter(codons)
#print (counts) #frequency of each codon


#read a table of codon score coefficient.
csc=open(sys.argv[2])
csc_=csc.readlines()
table={}
for line in csc_[1:]:
    line=line.strip()
    table[line.split("\t")[0]]=line.split("\t")[2]

print (table)
stop_codons=['TAA', 'TAG', 'TGA']
CSC_mean={}
for key in ref2:
    cscsum=0
    for codon in ref2[key]:
        if codon not in stop_codons:
            cscsum+=float(table[codon])
    CSC_mean[key]=cscsum/len(ref2[key])
print(list(CSC_mean.items())[1:5])

outFile=open(sys.argv[3],"w")
outFile.write("ID\tMean_CSC\n")
for key in CSC_mean:
    outFile.write(key+"\t"+str(CSC_mean[key])+"\n")
outFile.close()

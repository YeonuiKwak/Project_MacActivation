

''' Filter transcripts with 3UTR >1000 nt using Dic fuction.
Save last 500nt sequences
This could be applied to generate 3UTR sequence of specific region'''

'''usage:
python3 step6_6mercount_frontandlast500regions.py <custom3UTR.fa> <all6merIDandAUcontent.txt> <kmercountinfront500.txt> <kmercountinlast500.txt>
'''


import sys

inFile=open(sys.argv[1])
outFile_f=open(sys.argv[3],"w")
outFile_l=open(sys.argv[4],"w")
line = inFile.readlines()
dic={}
for i in range(0,len(line),2) :
    if line[i].strip()[0]==">":
        dic[line[i].strip()[1:]]=line[i+1]
#1.Front and Last 500 nt region sequence was extracted
front={}
last={}

for key, val in dic.items():
    if len(val)>=1000:
        #outFile.write(key+'\t'+val.strip()[-300:]+'\n')
        front[key]=val.strip()[:500]
        last[key]=val.strip()[-500:]

print(list(front.items())[1:5])


#2.

kmer_front_final={}
kmer_last_final={}



#count k-mers
for key in front:
    kmerdic={}
    i=0
    while i<(495):
        x=front[key][i:i+6]
        if x in kmerdic:
            kmerdic[x]+=1
        else:
            kmerdic[x]=1
        i=i+1

    kmer_front_final[key]=kmerdic

for key in last:
    kmerdic={}
    i=0
    while i<(495):
        x=last[key][i:i+6]
        if x in kmerdic:
            kmerdic[x]+=1
        else:
            kmerdic[x]=1
        i=i+1

    kmer_last_final[key]=kmerdic



#3. Load all 6-mer species
kmerlist=open(sys.argv[2])
kmer=[]
for x in kmerlist:
    kmer.append(x.split("\t")[0])
    #first column:6-mer species, secondcolumn : number of AU
    #count each of individaul 6-mer sequences in a given transcript.

outFile_f.write("id\t")
outFile_l.write("id\t")
for x in kmer:
    #Column name generation
    outFile_f.write(x+"\t")
    outFile_l.write(x+"\t")
outFile_l.write("\n")
outFile_f.write("\n")

#4. Build 6-mer count tables for the front 500 nucleotide regions.

#front region 500nt

for id in kmer_front_final:
    outFile_f.write(id+"\t")
    for a in kmer:
        if a not in kmer_front_final[id].keys():
            outFile_f.write(str(0)+"\t")
        else:
            outFile_f.write(str(kmer_front_final[id][a])+"\t")
    outFile_f.write("\n")

#last region 500nt

for id in kmer_last_final:
    outFile_l.write(id+"\t")
    for a in kmer:
        if a not in kmer_last_final[id].keys():
            outFile_l.write(str(0)+"\t")
        else:
            outFile_l.write(str(kmer_last_final[id][a])+"\t")
    outFile_l.write("\n")

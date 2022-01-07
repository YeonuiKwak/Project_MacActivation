


''' Filter transcripts with 3UTR >1000 nt using Dic fuction.
Save last 300nt sequences
This could be applied to generate 3UTR sequence of specific region'''

'''usage:
python3 UTR_countKmer.py ~/Desktop/TED-seq\ Rscript/3UTR_concat.fa ~/Desktop/ReadcountTable/Ref/all6merIDandAUcontent.txt kmercountinlast300.txt
'''


import sys

inFile=open(sys.argv[1])
outFile=open(sys.argv[3],"w")

line = inFile.readlines()
dic={}
for i in range(0,len(line),2) :
    if line[i].strip()[0]==">":
        dic[line[i].strip()[1:]]=line[i+1]
#300 nt regions.
front={}
last={}

for key, val in dic.items():
    if len(val)>=600:
        #outFile.write(key+'\t'+val.strip()[-300:]+'\n')
        front[key]=val.strip()[:300]
        last[key]=val.strip()[-300:]

print(list(front.items())[1:5])

kmer_front_final={}
kmer_last_final={}
#count k-mers
for key in front:
    kmerdic={}
    i=0
    while i<(295):
        x=front[key][i:i+6]
        if x in kmerdic:
            kmerdic[x]+=1
        else:
            kmerdic[x]=1
        i=i+1

    kmer_front_final[key]=kmerdic

for key in last:
    kmerdic={}
    while i<=(len(last[key])-6):
        if last[key][i:i+6] in kmerdic:
            kmerdic[last[key][i:i+6]]+=1
        else:
            kmerdic[last[key][i:i+6]]=0
        i=i+1

    kmer_last_final[key]=kmerdic


#print
#print(list(kmer_front_final.items())[1:5])

#kmerlist
kmerlist=open(sys.argv[2])
kmer=[]
for x in kmerlist:
    kmer.append(x.split("\t")[0])
#count each of individaul 6-mer sequences in a given transcript.
outFile.write("id\t")
for x in kmer:
    outFile.write(x+"\t")
outFile.write("\n")
'''
for id in kmer_front_final:
    outFile.write(id+"\t")
    for a in kmer:
        if a not in kmer_front_final[id]:
            outFile.write(str(0)+"\t")
        else:
            outFile.write(str(kmer_front_final[id][a])+"\t")
    outFile.write("\n")
'''


for id in kmer_last_final:
    outFile.write(id+"\t")
    for a in kmer:
        if a not in kmer_last_final[id]:
            outFile.write(str(0)+"\t")
        else:
            outFile.write(str(kmer_last_final[id][a])+"\t")
    outFile.write("\n")



outFile.close()

''' Filter transcripts with 3UTR >1000 nt using Dic fuction.
Save last 300nt sequences
This could be applied to generate 3UTR sequence of specific region'''


line = inFile.readlines()
dic={}
for i in range(0,len(line)) :
    l=line[i].strip()
    if l[0]==">":
        dic[l]=line[i+1]

for key, val in dic.items():
    if len(val)>=1000:
        outFile.write(key+'\n'+val.strip()[-300:]+'\n')

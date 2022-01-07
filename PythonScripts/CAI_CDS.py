

'''
flist=inFile.readlines()
count = len(flist)
x=0
while x < count:
    f=flist[x].strip()
    if f[0]==">":
        outFile.write(f+'\n')
        outFile.write(flist[x+1].strip())
        i=1
        while (x+2*i)<count:
                if flist[x].strip()==flist[x+2*i].strip():
                    outFile.write(flist[x+2*i+1].strip())
                    i=i+1
                else:
                    outFile.write("\n")
                    break
        x=x+2*i
    print x
outFile.close()
'''
'''
#Ensmeble  hg38 CDS sequences..
for line in inFile:
    line = line.strip()
    if line[0]== ">":
        l = line.split(' ')
        outFile.write('\n'+l[0]+'\n')
    else:
        outFile.write(line.strip())

outFile.close()
'''


'''
line = inFile.readlines()
dic={}
for i in range(0,len(line)) :
    l=line[i].strip()
    if l[0]==">":
        dic[l]=line[i+1].strip()

for key, val in dic.items():
    if  len(val)>302:
        outFile.write(key[1:]+'\t'+str(len(val))+'\n')
        outFile_second.write(key+'\n'+val+'\n')

outFile.close()
outFile_second.close()
'''

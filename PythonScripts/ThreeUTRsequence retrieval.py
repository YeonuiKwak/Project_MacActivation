

import sys
#Edit the fast sequence id.
'''
for example, #>ENST00000379389.4_utr3_1_0_chr1_1014479_f::chr1:1014478-1014540(+) will be converted to >ENST00000379389.4. 
'''
infile=open(sys.argv[1])
Outfile=open("3UTR.fa","w")
for line in infile:
    if line[0]==">":
        Id=line.strip().split("_")
        Outfile.write(Id[0]+"\n")
    else:
        Outfile.write(line)

Outfile.close()


#import sys
#Concatenate exon sequences under the same transcript id.
infile=open(sys.argv[1])
outfile=open("3UTR_concat.fa","w")

flist=infile.readlines()
count = len(flist)
x=0
while x < count:
    if flist[x][0]==">":
        outfile.write(flist[x])
        outfile.write(flist[x+1].strip())
        i=1
        while (x+2*i)<count:
                if flist[x]==flist[x+2*i]:
                    outfile.write(flist[x+2*i+1].strip())
                    i=i+1
                else:
                    outfile.write("\n")
                    break
        x=x+2*i
    print x
outfile.close()

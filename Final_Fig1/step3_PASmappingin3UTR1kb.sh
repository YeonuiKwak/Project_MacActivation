#cps.all.internalpAremoved.bed: only with positions of greater than 5 readcounts.
#usage: bash step3.sh <3UTR1kbext.bed> <CPS.all.internalpAremoved.gt5.bed>


bedtools intersect -b $1 -a $2 -s  -wao |\
awk 'BEGINS{OFS="\t"}$7!="."{print $1,$2,$3,$10":"$4,$5,$6}'> cpsin3UTR1kb.bed6
#assign CPS isoform id  to each CPS in a transcript.
awk '$6=="+"' cpsin3UTR1kb.bed6 | \
 sort -k4,4 -k2,2nr -u | \
 awk 'FNR==1{prev=$4;count=1;print $0"\t"count;next}
   {cur=$4;if(cur!=prev) count=1; else ++count; prev=cur; print $0"\t"count}' | \
 sort  -k1,1 -k4,4 -k2,2n -k3,3n -u > cpsin3UTR1kb.pl.bed7
awk '$6=="-"' cpsin3UTR1kb.bed6 | \
 sort -k4,4 -k2,2n -u | \
 awk 'FNR==1{prev=$4;count=1;print $0"\t"count;next}
   {cur=$4;if(cur!=prev) count=1; else ++count; prev=cur; print $0"\t"count}' | \
 sort -k1,1 -k4,4 -k2,2n -k3,3n -u > cpsin3UTR1kb.mn.bed7


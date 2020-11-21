#Step9: Assign transcript to CPS position.
#cps.all.internalpAremoved.bed: only with positions of greater than 5 readcounts.
bedtools intersect -b utr.1kbext.bed -a cps.all.internalpARemoved.bed -s -wao |\
awk 'BEGINS{OFS="\t"}$7!="."{print $1,$2,$3,$10,$5,$6}'> cpsin3UTR1kb.bed6
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
cat cpsin3UTR1kb.pl.bed7 cpsin3UTR1kb.mn.bed7 | \
	sort -k1,1 -k2,2n -k3,3n > cpsin3UTR1kb.bed7


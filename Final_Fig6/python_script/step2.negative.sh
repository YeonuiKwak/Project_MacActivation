#Negative regions: when distal UTR has shorter poly(A) tail.
#the region inclusive in this set of distal 3'UTRs shouldn't have the sequence motif found in the set of dsital 3' UTR with longer poly(A) tail .

bedtools subtract -a custom.3UTR.distal.shorter.bed12 -b custom.3UTR.proximal.longer.bed12 -s > Negative_distal.only.UTR.region.bed
#Extract The region in fastsa.
bedtools getfasta -fi ~/Work/shared/ref/hg38/hg38.fa -bed Negative_distal.only.UTR.region.bed -split -name -s | \
      fold -w 60 > Negative_custom3UTR.distalOnly.fa



#Extract custom 3'UTR sequence using custom.UTR.bed12.
#bedfile :contains 21224 UTRs: old version(before 1million CPS.bed normalization)
#bedtools getfasta -fi ~/Work/shared/ref/hg38/hg38.fa -bed custom.3UTR.bed12 -tab -split -s -name > custom3UTR.fa


#new custom 3UTR file:custom3UTR.n.30469.fa
#Extract custom 3'UTR sequence using custom.UTR.bed12.
#bedfile :contains 21224 UTRs: old version(before 1million CPS.bed normalization)
#bedtools getfasta -fi ~/Work/shared/ref/hg38/hg38.fa -bed custom.3UTR.n.30649.bed12 -tab -split -s -name > custom3UTR.n.30649.tab.fa

bedtools getfasta -fi ./hg38.fa -bed custom.3UTR.n.31490.bed12 -tab -split -s -name >custom3UTR.n.31490.tab.fa
bedtools getfasta -fi ./hg38.fa -bed custom.3UTR.n.31490.bed12 -split -s -name >custom3UTR.n.31490.fa

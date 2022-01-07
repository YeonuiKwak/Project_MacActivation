#kmer anlaysis
bedtools getfasta -fi ~/Work/shared/ref/hg38/hg38.fa -bed custom.3UTR.bed12 -name -split -s > custom3UTR.fa


#kmer analysis emploly \n separated fasta file.
~/Work/shared/code/find-motifs -i custom3UTR.n.30649.fa > customutr3.n.30649.kmer.txt

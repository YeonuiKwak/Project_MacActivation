#!/bin/bash
#use full window for each CPS position (the ouctome of merge of all CPS position with 5bp window), and
# calculate the coverage of each window by each CPS bam file.
#Do this for all biological replicates of time samples(n=8; 2replicates x 4timepoints )

#step1:


cat $1 $2 $3 $4\
|sort -k1,1 -k2,2n -k5,5nr\
|awk '$6=="+"{print $1"\t"$2-3"\t"$3+3"\t"$4"\t"$5"\t+";next}{print $1"\t"$2-3"\t"$3+3"\t"$4"\t"$5"\t-";}'\
|sort -k1,1 -k2,2n|bedtools merge -s -c 5,6 -o sum,distinct -i - \
|awk '$4>5{print $1"\t"$2"\t"$3"\t"$1":"int(($2+$3)/2)"\t"$4"\t"$5}' \
|sort -k1,1 -k2,2n >cps.all.full.window.internalpAremoved.gt5.bed

cat $1 $2 $3 $4\
|sort -k1,1 -k2,2n -k5,5nr\
|awk '$6=="+"{print $1"\t"$2-3"\t"$3+3"\t"$4"\t"$5"\t+";next}{print $1"\t"$2-3"\t"$3+3"\t"$4"\t"$5"\t-";}'\
|sort -k1,1 -k2,2n|bedtools merge -s -c 5,6 -o sum,distinct -i - \
|awk '$4>5{print $1"\t"int(($2+$3)/2)"\t"int(($2+$3)/2)+1"\t"$1":"int(($2+$3)/2)"\t"$4"\t"$5}' \
|sort -k1,1 -k2,2n >cps.all.PASpositions.internalpAremoved.gt5.bed

#!/bin/bash
#step4. Convert fitlered.fasta id to bedfile
sample=(
0h
1h
2h
4h
)
for i in ${sample[@]}; do
cat cps.${i}.10ntwindow_filtered.pl.fa.fasta|awk '/^>/{a=substr($1,2);print a}'|\
cut -d':' --output-delimiter=" " -f1,2 |\
cut -d'-' --output-delimiter=" " -f1,2 |\
awk 'BEGIN{OFS="\t"}{print $1,$3-1,$3,$1":"$3-1,".","+"}'>cps.${i}.filtered.pl.bed
cat cps.${i}.10ntwindow_filtered.mn.fa.fasta|\
awk '/^>/{a=substr($1,2);print a}'|\
cut -d':' --output-delimiter=" " -f1,2 |\
cut -d'-' --output-delimiter=" " -f1,2 |\
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$1":"$2-1,".","-"}'>cps.${i}.filtered.mn.bed 
cat cps.${i}.filtered.pl.bed cps.${i}.filtered.mn.bed >cps.${i}.filtered.bed
done

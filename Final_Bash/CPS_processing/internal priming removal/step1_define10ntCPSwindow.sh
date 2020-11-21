#!/bin/bash
#step1. make +/- 10 nt window of the CPS position in a bed file.
sample=(
0h
1h
2h
4h
)
for i in ${sample[@]}; do
awk '$6=="+"{print $1"\t"$2-9"\t"$3"\t"$4"\t"$5"\t"$6}' cps.${i}.bed>cps.${i}.10ntwindow.pl.bed
awk '$6=="-"{print $1"\t"$2"\t"$3+9"\t"$4"\t"$5"\t"$6}' cps.${i}.bed>cps.${i}.10ntwindow.mn.bed
done

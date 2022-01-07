#!/bin/bash

SAMPLE=(
0h
1h
2h
4h
)



for s in ${SAMPLE[@]}; do
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' cps.${s}.withoutinternalpA.1million.bed|sort -k1,1 -k2,2n >tmp.bed &
wait
bedtools map -a $1 -b tmp.bed -s -o sum -c 5 | cut -f1-7 |awk '$7>0{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
|sort -k1,1 -k2,2n >cps.rc.${s}.txt
done

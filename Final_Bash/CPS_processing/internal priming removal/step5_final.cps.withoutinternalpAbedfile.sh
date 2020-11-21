#!/bin/bash
#step5. Intersect with cps.0h.bed vs cps.0h.filtered.bed
sample=(
0h
1h
2h
4h
)
for i in ${sample[@]}; do
  bedtools intersect -a ~/Work2/yk724/CPSseq/CPS_positions/bed/cps.${i}.bed -b cps.${i}.filtered.bed -s  -wao |\
  awk 'BEGINS{OFS="\t"}$7!="."{print $1,$2,$3,$4,$5,$6}'|sort -k1,1 -k2,2 -k5,5nr -k6,6 -u > cps.${i}.withoutinternalpA.bed
done

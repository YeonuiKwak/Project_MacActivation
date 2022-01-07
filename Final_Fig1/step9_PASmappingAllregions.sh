#!/bin/bash
#usage: $1: transcript.bed13 $2: full.window.gt5,bed
cat $1|awk '$6=="+"{print $1,$2,$3+1000,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13;next}{if ($2>1000) print $1,$2-1000,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13;else print $0;next}'|
bedtools intersect -b $1 -a $2 -s -wao |\
awk 'BEGINS{OFS="\t"}$7!="."{print $1,$2,$3,$10":"$4,$5,$6,$19}'>cpsinallregions.bed

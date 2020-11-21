#!/bin/bash
#retrieve sequence of 10 nt window of CPS positions from the hg38.fa 
sample=(
0h
1h
2h
4h
)

#step2. retrieve sequences within the window&  convert lower case to uppercase.
for i in ${sample[@]}; do
seqtk subseq ~/Work/shared/ref/hg38/hg38.fa cps.${i}.10ntwindow.pl.bed >cps.${i}.10ntwindow.pl.fa
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' cps.${i}.10ntwindow.pl.fa > tmp
mv tmp cps.${i}.10ntwindow.pl.fa
seqtk subseq ~/Work/shared/ref/hg38/hg38.fa cps.${i}.10ntwindow.mn.bed >cps.${i}.10ntwindow.mn.fa
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' cps.${i}.10ntwindow.mn.fa > tmp
mv tmp cps.${i}.10ntwindow.mn.fa
done


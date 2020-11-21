#!/bin/bash

sample=(
0h
1h
2h
4h
)
#step3 Filter out the regions that have cosecutive AAAAA.
echo "filterout statistics">./prinseqinternalprimingfilter.gd 
for i in ${sample[@]}; do
prinseq-lite.pl -fasta cps.${i}.10ntwindow.pl.fa -out_format 1 -custom_params "A 5" -out_good cps.${i}.10ntwindow_filtered.pl.fa -out_bad cps.${i}.bad.pl.fa 2>> ./prinseqinternalprimingfilter.gd
prinseq-lite.pl -fasta cps.${i}.10ntwindow.mn.fa -out_format 1 -custom_params "T 5" -out_good cps.${i}.10ntwindow_filtered.mn.fa -out_bad cps.${i}.bad.mn.fa 2>> ./prinseqinternalprimingfilter.gd
done


#!/bin/bash
sample=(
  0h
  1h
  2h
  4h
  )

for i in ${sample[@]}; do
awk '$6=="+"{print $1"\t"$2"\t"$3"\t"$5}' cps.${i}.withoutinternalpA.bed>cps.${i}.withoutinternalpA.pl.bedGraph 
awk '$6=="-"{print $1"\t"$2"\t"$3"\t"$5}' cps.${i}.withoutinternalpA.bed>cps.${i}.withoutinternalpA.mn.bedGraph 
done

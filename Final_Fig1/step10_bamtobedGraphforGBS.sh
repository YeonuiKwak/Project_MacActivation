#step10: For Genome browser shot, Convert bed file to bam file
#!/bin/bash
sample=(
  0h
  1h
  2h
  4h
  )

for i in ${sample[@]}; do
awk '$6=="+"{print $1"\t"$2"\t"$3"\t"$5}' cps.${i}.withoutinternalpA.1million.bed>cps.${i}.withoutinternalpA.1million.pl.bedGraph 
awk '$6=="-"{print $1"\t"$2"\t"$3"\t"$5*(-1)}' cps.${i}.withoutinternalpA.1million.bed>cps.${i}.withoutinternalpA.1million.mn.bedGraph 
done

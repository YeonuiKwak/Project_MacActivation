cat hg38_CDS.fa |paste - -|awk 'BEGIN{FS="_";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7}'|sort -k1,1 -k3,3n - |awk 'NR==1{prev=$1;count=1;seq[prev]=$8;next}{cur=$1;if(cur==prev) {++count;seq[prev]=seq[prev]""$8;next}else{print prev"\t"seq[prev]"\t"count;count=1;prev=cur;seq[prev]=$8;next}}'>hg38_CDS.concat.fa


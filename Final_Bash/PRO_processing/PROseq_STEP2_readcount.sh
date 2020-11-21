# Gene body counts
awk -v range=500 '($3-$2)>1000{print $1"\t"$2+range"\t"$3-range"\t"$4"\t"$5"\t"$6;next}' $1 | sort -k1,1 -k2,2n -k3,3n > _gb.bed.tmp
bedtools map -a _gb.bed.tmp -b $2 -o sum  -c 4 | cut -f4,6,7 > _pl.gb.tmp
bedtools map -a _gb.bed.tmp -b $3 -o sum -c 4 | cut -f4,7 > _mn.gb.tmp
# Promoter proximal counts
awk '$2>100{print $0}' $1|awk -v range=500 '$6=="+"{print $1"\t"$2-100"\t"$2+range-100"\t"$4"\t"$5"\t"$6;next}$3>range{print $1"\t"$3-range+100"\t"$3+100"\t"$4"\t"$5"\t"$6}'| sort -k1,1 -k2,2n > _pp.bed.tmp
bedtools map -a _pp.bed.tmp -b $2 -o sum -c 4 | cut -f4,6,7 > _pl.pp.tmp
bedtools map -a _pp.bed.tmp -b $3 -o sum -c 4 | cut -f4,7 > _mn.pp.tmp

paste _pl.gb.tmp _mn.gb.tmp | awk '$2=="+"{print $1"\t"0+$3;next}{print $1"\t"0+$5}' > _gb.tmp
paste _pl.pp.tmp _mn.pp.tmp | awk '$2=="+"{print $1"\t"0+$3;next}{print $1"\t"0+$5}' > _pp.tmp

awk 'BEGIN{print "id\tgb\tpp"}NR==FNR{pp[$1]=$2;next}{print $1"\t"$2"\t"0+pp[$1]}' _pp.tmp _gb.tmp
# rm _??.*.tmp

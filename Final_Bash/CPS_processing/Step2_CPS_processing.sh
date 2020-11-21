#Step1: Pre-processing CPS-seq: Fastq to Bed and BedGraph

sample=(








)
for i in ${sample[@]}; do
bash ~/Work2/yk724/TEDseq/CPS/bash/Step0_CPSseq_preprocessing.sh -genome ~/Work/shared/ref/hg38/ -c ~/Work/shared/ref/hg38/chrNameLength.txt -fastq fastq/${i}
done



#step2: Removal internal priming candidates

2.0: Collapse replicate 1 and replicate 2 samples

#<step1_collapse_CPS.sh> Code below:
#!/bin/bash
#usage: bash step1_collapse_CPS.sh <biological rep1.bed> <biological rep2.bed> <outputName.bed>
#collapase cps positions betwee two biological replicates of CPS-seq.

cat $1 $2 \
|sort -k1,1 -k2,2n -k5,5nr |sort -k1,1 -k2,2n -u > $3
######

#Execulte the script below.
bash ../bash/step1_collapse_CPS.sh 0h.r1.bed6 0h.r2.bed6 0h.all.bed6
bash ../bash/step1_collapse_CPS.sh 1h.r1.bed6 1h.r2.bed6 1h.all.bed6
bash ../bash/step1_collapse_CPS.sh 2h.r1.bed6 2h.r2.bed6 2h.all.bed6
bash ../bash/step1_collapse_CPS.sh 4h.r1.bed6 4h.r2.bed6 4h.all.bed6



2.1: Write 10 nucleotide windows "downstream " of CPS position.

sample=(
0h
1h
2h
4h
)
for i in ${sample[@]}; do
  awk '$6=="+"{print $1"\t"$2"\t"$2+10"\t"$4"\t"$5"\t"$6}' cps.${i}.bed>cps.${i}.10ntwindow.pl.bed
  awk '$6=="-"{print $1"\t"$3-10"\t"$3"\t"$4"\t"$5"\t"$6}' cps.${i}.bed>cps.${i}.10ntwindow.mn.bed
done


2.2: Retrieve sequece from 10ntwindow.pl/mn.bed

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
bedtools getfasta ~/Work/shared/ref/hg38/hg38.fa -bed cps.${i}.10ntwindow.pl.bed -split -name -s >cps.${i}.10ntwindow.pl.fa
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' cps.${i}.10ntwindow.pl.fa > tmp
mv tmp cps.${i}.10ntwindow.pl.fa

bedtools getfasta~/Work/shared/ref/hg38/hg38.fa -bed cps.${i}.10ntwindow.mn.bed -split -name -s >cps.${i}.10ntwindow.mn.fa
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' cps.${i}.10ntwindow.mn.fa > tmp
mv tmp cps.${i}.10ntwindow.mn.fa
done


#step3. Filterout the regions that have consecutive As.

echo "filterout statistics">./prinseqinternalprimingfilter.gd

for i in ${sample[@]}; do
prinseq-lite.pl -fasta cps.${i}.10ntwindow.pl.fa -out_format 1 -custom_params "A 5" -out_good cps.${i}.10ntwindow_filtered.pl -out_bad cps.${i}.bad.pl 2> ./prinseqinternalprimingfilter.gd
prinseq-lite.pl -fasta cps.${i}.10ntwindow.mn.fa -out_format 1 -custom_params "A 5" -out_good cps.${i}.10ntwindow_filtered.mn -out_bad cps.${i}.bad.mn 2>> ./prinseqinternalprimingfilter.gd
done

#!/bin/bash
#step4. Convert filtered.fasta id to bedfile
#reference: http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
#note:fasta id is start site has +1 format, thus start site should be shifted by -1 back to normal bed format.
sample=(
0h
1h
2h
4h
)
for i in ${sample[@]}; do
cat cps.${i}.10ntwindow_filtered.pl.fa.fasta|awk '/^>/{a=substr($1,2);print a}'|cut -d':' --output-delimiter=" " -f1,2 |cut -d'-' --output-delimiter=" " -f1,2|awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,$1":"$2-1,".","+"}'>cps.${i}.filtered.pl.bed
cat cps.${i}.10ntwindow_filtered.mn.fa.fasta|awk '/^>/{a=substr($1,2);print a}'|cut -d':' --output-delimiter=" " -f1,2 |cut -d'-' --output-delimiter=" " -f1,2|awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,$1":"$3-1,".","-"}'>cps.${i}.filtered.mn.bed
cat cps.${i}.filtered.pl.bed cps.${i}.filtered.mn.bed >cps.${i}.filtered.bed
done



#step5. Intersect with cps.0h.bed vs cps.0h.filtered.bed
#Add the readcount information.
sample=(
0h
1h
2h
4h
)
for i in ${sample[@]}; do
  bedtools intersect -a cps.${i}.bed -b cps.${i}.filtered.bed -s  -wao |\
  awk 'BEGINS{OFS="\t"}$7!="."{print $1,$2,$3,$4,$5,$6}'> cps.${i}.withoutinternalpA.bed
done

#Before step6. Normalize the read count to 1 million per time point.
#Step6:Collapse all temporal CPS.beds to all.CPS.bed (all.CPS.internalpA.removed.bed)
#Step6.1:set the full CPS window, and sum readcounts mapped in the windows, and define the midpoint as PAS.
#input:cps.0h|1h|2h|4h.withoutinternlapA.bed (normalized to 1 million)
#output: "cps.all.internalpAremoved.bed"
#output: "cps.all.full.window.internalpAremoved.bed"


cat $1 $2 $3 $4\
|sort -k1,1 -k2,2n -k5,5nr\
|awk '$6=="+"{print $1"\t"$2-3"\t"$3+3"\t"$4"\t"$5"\t+";next}{print $1"\t"$2-3"\t"$3+3"\t"$4"\t"$5"\t-";}'\
|sort -k1,1 -k2,2n|bedtools merge -s -c 5,6 -o sum,distinct -i - \
|awk '$4>0{print $1"\t"$2"\t"$3"\t"$1":"int(($2+$3)/2)"\t"$4"\t"$5}' \
|sort -k1,1 -k2,2n >cps.all.full.window.internalpAremoved.bed

#Step7:In R, filter out CPSs with less than 5 readcounts.->Save the file





#Step8:In R, Extend last exon of 3'UTR by 1kb.->"utr.1kbext.bed"
#Step9: Assign transcript to CPS position.->"cpsin3UTR1kb.bed7"
#Input for step9 : cps.all.internalpAremoved.bed: only with positions of greater than 5 readcounts.
bedtools intersect -b utr.1kbext.bed -a cps.all.internalpARemoved.bed -s  -wao |\
awk 'BEGINS{OFS="\t"}$7!="."{print $1,$2,$3,$10,$5,$6}'> cpsin3UTR1kb.bed6
#assign CPS isoform id  to each CPS in a transcript.
awk '$6=="+"' cpsin3UTR1kb.bed6 | \
	sort -k4,4 -k2,2nr -u | \
	awk 'FNR==1{prev=$4;count=1;print $0"\t"count;next}
		{cur=$4;if(cur!=prev) count=1; else ++count; prev=cur; print $0"\t"count}' | \
	sort  -k1,1 -k4,4 -k2,2n -k3,3n -u > cpsin3UTR1kb.pl.bed7
awk '$6=="-"' cpsin3UTR1kb.bed6 | \
	sort -k4,4 -k2,2n -u | \
	awk 'FNR==1{prev=$4;count=1;print $0"\t"count;next}
		{cur=$4;if(cur!=prev) count=1; else ++count; prev=cur; print $0"\t"count}' | \
	sort -k1,1 -k4,4 -k2,2n -k3,3n -u > cpsin3UTR1kb.mn.bed7
cat cpsin3UTR1kb.pl.bed7 cpsin3UTR1kb.mn.bed7 | \
	sort -k1,1 -k2,2n -k3,3n > cpsin3UTR1kb.bed7



#Step10: After step5, Normalize <cps.${i}.withoutinternalpA.bed> the total read count per sample to 1 million in R. Make strand-specific bedGraph file. For Genomebrowser shot.

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



#Step11: Scatterplot between two biological replicates.
Repeat step5, but use cps.0h.rep1|rep2.bed instead of merged file <cps.0h.bed>
sample=(
0h.rep1
0h.rep2
1h.rep1
1h.rep2
2h.rep1
2h.rep2
4h.rep1
4h.rep2
)
for i in ${sample[@]}; do
  bedtools intersect -a cps.${i}.bed -b cps.${i}.filtered.bed -s  -wao |\
  awk 'BEGINS{OFS="\t"}$7!="."{print $1,$2,$3,$4,$5,$6}'> cps.${i}.withoutinternalpA.bed
done


#Step12: Finally, Per each PAS window, calculate sum of readcounts in each sample. (n=8)

#code usage:
bash ../step2.CalculateReadcount.sh cps.all.full.window.internalpAremoved.bed
input: cps.all.full.window.internalpAremoved.gt5.bed or cps.all.full.window.internalpAremoved.bed
#build readcount table for the regions without internal priming.
#!/bin/bash
#After filtering cps.position window with readcounts >5.
#2)Calculate the coverage of this region by CPS-peak reads.
#Do this for all biological replicates of time samples(n=8; 2replicates x 4timepoints )


SAMPLE=(
0h.r1
0h.r2
1h.r1
1h.r2
2h.r1
2h.r2
4h.r1
4h.r2
)




for s in ${SAMPLE[@]}; do
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' cps.${s}.withoutinternalpA.bed|sort -k1,1 -k2,2n >tmp.bed &
wait
bedtools intersect -wao -s -a $1 -b tmp.bed \
|awk '$7!="."{print $0}'|cut -f1-6,11 -|bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,sum -i - \
|awk '$7>0{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
|sort -k1,1 -k2,2n >cps.rc.${s}.txt

done

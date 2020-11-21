
#PROseq converted <bam> to <bed>
bedtools bamtobed -i $1 | awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' |      awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,"-"}; ($6 == "-") {print $1,$3-1,$3,$4,$5,"+"}' | gzip >$2.bed.gz

#remove reads mapped to mitochondria chromosomes, and rRNA"
gunzip -c $2.bed.gz |grep "rRNA|chrM" -v |grep  "_" -v|sort-bed - |gzip >$2.nr.rs.bed.gz

#sort by chromosome number and location within the chromosome.
gunzip -c $2.nr.rs.bed.gz |sort -k1,1 -k2,2n -k3,3n> sorted.$2.nr.rs.bed
#Convert bed to bedGraph.
bedtools genomecov -bg -i sorted.$2.nr.rs.bed -g ~/Work/shared/ref/hg38/chrNameLength.txt -strand + >$2.pl.bedGraph
bedtools genomecov -bg -i sorted.$2.nr.rs.bed -g ~/Work/shared/ref/hg38/chrNameLength.txt -strand - >$2.mn.bedGraph
head $2.bedGraph

#../../scripts/proseq_coverage.sh ../../../RefAnno/hg19GencodeAnnotation_genes.bed THP1.rep2.d.0.pl.bedGraph THP1.rep2.d.0.mn.bedGraph >THP1.rep2.d.0.txt
#../../scripts/proseq_coverage.sh ../../../RefAnno/sorted_hg19GencodAnno.bed THP1.rep2.d.0.pl.bedGraph THP1.rep2.d.0.mn.bedGraph >THP1.rep2.d.0.txt
#head THP1.rep2.d.0.txt

#Check the reads for TNF.
#cat THP1.rep2.d.0.txt|grep "ENST0000044926"

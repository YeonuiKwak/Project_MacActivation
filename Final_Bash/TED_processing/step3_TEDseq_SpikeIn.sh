#Fastq file pre-processing
#set working directory to yk724/TEDseq/
## Spike-in RNA alignment and check the length distribution of spike-in RNAs.
while test $# -gt 0; do
	case "$1" in
		-h|--help)
	echo ""
	echo "Post-QC TEDseq.fastq will be aligned to spike-in genome"
	echo ""
	echo "Take PREFIX.QC_end.fastq.gz"
	echo "or  *.QC_end.fastq can be used as an input."
	echo "Output format is  sort.bam, and bedGraph/bigWig."
	echo "Usage: bash TEDseq_SpikeIn.sh <InputFastqFile>"
			exit 0
			;;
		-I)
			shift
			if test $# -gt 0; then
				export FQ_INPUT=$1
			else
				echo "no FASTQ input is specified"
				exit 1
			fi
			shift
			;;
	esac
done
echo ""
echo "Post-QC:TED-seq reads ->Spike-in alignment"
#bwa index ~/Work/labuser/yk724/TEDseq/spike-in/spikein.fa
export name=`echo ${FQ_INPUT}|rev|cut -d \/ -f 1|cut -d\. -f 3-|rev`
export TMP="./${name}_SpikeIn"
export BWAIDX=~/Work/labuser/yk724/TEDseq/spike-in/spikein.fa
export CHINFO=~/Work/labuser/yk724/TEDseq/spike-in/spike_Len.txt
mkdir ${TMP}
mkdir ${TMP}/out
echo "QC_end.fastq; number of reads:" > ${TMP}/Align.log
gunzip -c ${FQ_INPUT}|grep @ -c >>${TMP}/Align.log
gunzip -c ${FQ_INPUT}>${TMP}/${name}.fastq &
wait
echo "START bwa -mem alignment." 
bwa mem -t 8 ~/Work/labuser/yk724/TEDseq/spike-in/spikein.fa ${TMP}/${name}.fastq -f ${TMP}/out/${name}_spikeIn.sam &
wait
echo ""
echo "START convert sam to bam. sort bam."

samtools view  -b -q 0 ${TMP}/out/${name}_spikeIn.sam| \
samtools sort -n -@ 8 - >${TMP}/out/${name}_spikeIn_sort.bam
#bedfiles of mapping 5' end of reads to spike-in
bedtools bamtobed -i ${TMP}/out/${name}_spikeIn_sort.bam 2> kill.warnings|awk 'BEGIN{OFS="\t"};\
($5>0){print $0}' |awk 'BEGIN{OFS="\t"};($6=="+"){print $1,$2,$2+1,$4,$5,$6};\
($6=="-"){print $1,$3-1,$3,$4,$5,$6}' |gzip >${TMP}/out/${name}_spikein.bed.gz
echo "spike-In  mapping  readcount and  positions" 
echo "Spike-in : readcounts per poly(A) length standards" >>${TMP}/Align.log
echo "total mappable readcount:"
gunzip -c ${TMP}/out/${name}_spikein.bed.gz |wc -l >>${TMP}/Align.log
echo "mappable reads per poly(A) length standard:"
gunzip -c ${TMP}/out/${name}_spikein.bed.gz |cut -f 1|sort|uniq -c >>${TMP}/Align.log
echo "A40 length:">>${TMP}/Align.log
gunzip -c ${TMP}/out/${name}_spikein.bed.gz |awk '$1=="A40"{a +=$2;b +=1}END{print a/b}'>>${TMP}/Align.log
echo "A80 length:">>${TMP}/Align.log
gunzip -c ${TMP}/out/${name}_spikein.bed.gz |awk '$1=="A80"{a +=$2;b +=1}END{print a/b}'>>${TMP}/Align.log
echo "A120 length:">>${TMP}/Align.log
gunzip -c ${TMP}/out/${name}_spikein.bed.gz |awk '$1=="A120"{a +=$2;b +=1}END{print a/b}'>>${TMP}/Align.log
echo "A160 length:">>${TMP}/Align.log
gunzip -c ${TMP}/out/${name}_spikein.bed.gz |awk '$1=="A160"{a +=$2;b +=1}END{print a/b}'>>${TMP}/Align.log
#bigWig
## Convert to bedGraph ... Can't gzip these, unfortunately.
#bed should be sorted before  creating bedGraph.
gunzip -c ${TMP}/out/${name}_spikein.bed.gz|sort-bed -|gzip >${TMP}/out/${name}_spikein.sorted.bed.gz
bedtools genomecov -bg -i ${TMP}/out/${name}_spikein.sorted.bed.gz -g ~/Work/labuser/yk724/TEDseq/spike-in/spike_Len.txt -strand + > ${TMP}/out/${name}.plus.bedGraph


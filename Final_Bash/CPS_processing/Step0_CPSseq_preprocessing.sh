 #!/usr/bin/env bash
while test $# -gt 0
do
        case "$1" in
            -h|-help)
                echo ""
                echo "Preprocesses CPS_annotation RNA-seq data."
                echo "usage: CPS_annotation.sh -genome ~/Work/shared/ref/hg38/ -c ~/Work/shared/ref/hg38/chrNameLength.txt -fastq ./*.fastq.gz"
                echo "Prefix *(.fastq.gz) should be used."
                echo "or *.fastq.gz in the current working directory as input and writes"
                echo "BAM and bigWig files as output to the user_Assigned output directory"
                echo "requirement in the current working directory: cutadapt 1.8.3, fastx_trimmer, seqtk, prinseq-lite.pl 0.20.2, bwa, samtools, bedtools, and bedGraphToBigWig."
                echo "set PATH =$PATH:/home/miniconda/bin/"
                echo "--UMI1=8 [The length of UMI barcode on the 5' end of R1 Read]"
                exit 0
                ;;
            -genome)
                shift
                if test $# -gt 0; then
                    export STARIDX=$1
                else
                    echo "no STAR Index Directory specified"
                    exit 1
                fi
                shift
                ;;
            -c)
                shift
                if test $# -gt 0; then
                    export CHINFO=$1
                else
                    echo "CHINFO is not specified"
                    exit 1
                fi
                shift
                ;;
            -fastq)
                shift
                if test $# -gt 0; then
                    export  FQ_INPUT=$1
                else
                    echo "FQ_INPUT is not specified"
                    exit 1
                fi
                shift
                ;;

        esac
done


#Check arguments
if [ -z "$FQ_INPUT" ]; then
    echo "-fastq is required."
    exit 1
fi
if [ -z "$STARIDX" ]; then
    echo "star-index is required."
    exit 1
fi
if [ -z "$CHINFO" ]; then
    echo "-c is required"
    exit 1
fi




#Preprocess CPS_annotation RNA-seq data
if [ -z "${UMI1}" ]; then
	export UMI1=8
fi
if [ -z "${UMI2}" ];then
    export UMI2=0
fi
echo " "
echo "processing CPS annotation Seq data ..."
echo " "
echo "SEQ         CPS Annotation"


echo "Report 5' ends"
echo "reporting 5' end of the reads in the opposite strand is mapping 3' cleavage site of the corresponding RNAs."
echo ""
echo "input/file paths:"
echo "star indexed genome directory         $STARIDX"
echo "chromInfo                            $CHINFO"

#Make directories
TMPDIR="./tmp"
mkdir ${TMPDIR}
mkdir ${TMPDIR}/noadapt
mkdir ${TMPDIR}/passQC
dedup_l=30
UMI1=8
OligodT=8
n1=$[UMI1+OligodT]
mkdir ${TMPDIR}/noadapt/l${dedup_l}
mkdir ${TMPDIR}/noadapt/l${dedup_l}_nodups
Name=`echo ${FQ_INPUT}| rev| cut -d \/ -f 1| rev| cut -d \. -f 1`
OUTPUT="./${Name}_output"
mkdir ${OUTPUT}

echo ""
echo "preprocessing fastq files"
##read stats
echo "${Name}.fastq" >${TMPDIR}/${Name}.QC.log
echo "Number of original input reads:" >>${TMPDIR}/${Name}.QC.log
gunzip -c ${FQ_INPUT}.fastq.gz| grep @ -c >>${TMPDIR}/${Name}.QC.log &
wait
###############################################################################
#gunzip the fastq.gz file and save it in the current directory.
gunzip -c ${FQ_INPUT}.fastq.gz >${Name}.fastq &
wait

echo "adaptor dimer removal"
if [ -z "$thread" ]; then
   thread=8
fi

if [ -z "$ADAPT_SE" ]; then
  ADAPT_SE="TGGAATTCTCGGGTGCCAAGG" #equal to RA3 sequence
fi

## Remove adapter, UMI barcode, additional barcode, and low quality (q=20) base from 3prime end of reads. Keep read length >=15 after trimmming
# Remove adapter
  cutadapt -a ${ADAPT_SE} -e 0.10 --overlap 5 --minimum-length=30 --nextseq-trim 20 --output=${TMPDIR}/${Name}_q20trim.fastq --untrimmed-output=${TMPDIR}/${Name}_q20untrim.fastq ${FQ_INPUT}.fastq.gz 2>${Name}_cutadaptreport.txt

wait


#By using the --cut option or its abbreviation -u, it is possible to unconditionally remove bases from the beginning or end of each read. If the given length is positive, the bases are removed from the beginning of each read. If it is negative, the bases are removed from the end.
#3' end of Read1  doesn't have any other additional sequence than adaptor.




#minimum length =30, phred quality score>=20 for both trimmed and untrimmed.
#cutadapt --minimum-length=30 ${TMPDIR}/${Name}_untrim.fastq --output=${TMPDIR}/${Name}_q20untrim.fastq --nextseq-trim 20 &
#wait

#cutadapt --minimum-length=30 ${TMPDIR}/${Name}_trim.fastq --output=${TMPDIR}/${Name}_q20trim.fastq --nextseq-trim 20 &
#wait

#sort by fastq id.
cat ${TMPDIR}/${Name}_q20untrim.fastq ${TMPDIR}/${Name}_q20trim.fastq | paste - - - - |LC_ALL=C sort --temporary-directory=${TMPDIR} --parallel=10 -k1,1 -S 10G | tr '\t' '\n' > ${TMPDIR}/noadapt/${Name}_noadapt.fastq &
wait


#
echo 'Number of reads after adapter removal and QC:' >> ${TMPDIR}/${Name}.QC.log
cat ${TMPDIR}/noadapt/${Name}_noadapt.fastq | grep @ -c >> ${TMPDIR}/${Name}.QC.log

#Remove the old fastq. we only need ${TMPDIR}/noadapt/${Name}_noadapt.fastq
#rm ${TMPDIR}/${Name}_trim.fastq ${TMPDIR}/${Name}_untrim.fastq ${TMPDIR}/${Name}_q20trim.fastq




  if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then
   # Remove PCR duplciates.
   # fastq with the first dedup_l nt
   cat ${TMPDIR}/noadapt/${Name}_noadapt.fastq | fastx_trimmer -Q33 -l ${dedup_l} -o ${TMPDIR}/noadapt/l${dedup_l}/${Name}_noadapt_l${dedup_l}.fastq &
   wait
   # deduplicate using the first dedup_L nt
   prinseq-lite -fastq ${TMPDIR}/noadapt/l${dedup_l}/${Name}_noadapt_l${dedup_l}.fastq -derep 1 -out_format 3 -out_bad null -out_good ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup -min_len 8 2> ${OUTPUT}/${Name}.prinseq-pcrDups.gd &
   wait
   rm ${TMPDIR}/noadapt/l${dedup_l}/${Name}_noadapt_l${dedup_l}.fastq &
   wait

   # make a list of name from deduplicated fastq
   cat ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup.fastq | awk '(NR%4==1){print substr($1,2)}' > ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}.txt
   wait

   rm ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup.fastq

   # generate fastq from the list of name
   seqtk subseq ${TMPDIR}/noadapt/${Name}_noadapt.fastq ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}.txt > ${TMPDIR}/passQC/${Name}_dedup_withBarcode.fastq &
   wait
   rm ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}.txt

   # trim the UMI and additional barcode after dereplicate
   #8 T and 8N UMI.
   prinseq-lite -trim_left ${n1} -fastq ${TMPDIR}/passQC/${Name}_dedup_withBarcode.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${Name}_dedup_BarcodeRemoved 2>> ${OUTPUT}/${Name}.prinseq-pcrDups.gd &
   wait
   # min_len 15
   prinseq-lite -min_len 15 -fastq ${TMPDIR}/passQC/${Name}_dedup_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${Name}_dedup_QC_end 2>> ${OUTPUT}/${Name}.prinseq-pcrDups.gd &
   wait

   echo 'Number of reads after PCR duplicates removal and QC:' >> ${TMPDIR}/${Name}.QC.log
   cat ${TMPDIR}/passQC/${Name}_dedup_QC_end.fastq | grep @ -c >> ${TMPDIR}/${Name}.QC.log

  fi

  cat ${OUTPUT}/${Name}.prinseq-pcrDups.gd

wait
rm -r ${TMPDIR}/noadapt
rm ${TMPDIR}/passQC/${Name}_dedup_withBarcode.fastq ${TMPDIR}/passQC/${Name}_dedup_BarcodeRemoved.fastq &
wait
#gzip ${TMPDIR}/passQC/${Name}_dedup_QC_end.fastq

#########################Alignment to genome###########################

echo " "
echo "Mapping reads:"

#QC_INPUT=`ls  ${TMPDIR}/passQC/*_QC_end.fastq | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- | cut -d _ -f 2- | rev| sort | uniq`

#align using star
~/Work/bin/star --genomeDir $STARIDX \
        --readFilesIn ${TMPDIR}/passQC/${Name}_dedup_QC_end.fastq \
        --outFilterMultimapNmax 1 --runThreadN 8 \
        --outFileNamePrefix ${TMPDIR}/passQC/${Name}. &

wait
samtools view -b -q 10 ${TMPDIR}/passQC/${Name}.Aligned.out.sam|\
samtools sort -@ 8 - >${TMPDIR}/${Name}.sort.bam &
wait
########################
rm ${TMPDIR}/passQC/${Name}.Aligned.out.sam
cp ${TMPDIR}/${Name}.sort.bam ./${OUTPUT} & ## Saves the sorted BAM in the output file.
wait


## Cleanup
#Remove bam file in the temporary directory if their size is less than 1024 kilobytes.
find ${TMPDIR} -name "*.sort.bam" -size -1024k -delete


#############################################
## Write out the bigWigs.
echo " "
echo "Writing bigWigs:"

    f=`echo "${TMPDIR}/${Name}.sort.bam"`
    j=`echo $f | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev `
    echo $j > ${OUTPUT}/${j}.align.log

# in SE, mapping 5' end position of the reads on the opposite strand.


#if [[ "${RNA5}" == "R1_5prime" && "${OPP}" == "FALSE" ]] ; then ## report The 5 prime end of the RNA.   #like GRO-seq
#if [[ "$SE_OUTPUT" == "G" ]] ; then
#  bedtools bamtobed -i $f 2> ${TMPDIR}/kill.warnings| awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \

# awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | gzip > ${TMPDIR}/$j.bed.gz
#elif [[ "${RNA3}" == "R1_5prime" && "${OPP}" == "TRUE" ]] ; then  #like PRO-seq

bedtools bamtobed -i $f 2> ${TMPDIR}/kill.warnings| awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,"-"}; ($6 == "-") {print $1,$3-1,$3,$4,$5,"+"}' | gzip > ${TMPDIR}/$j.bed.gz


echo 'Number of mappable reads:' >> ${OUTPUT}/${j}.align.log
   echo `gunzip -c ${TMPDIR}/$j.bed.gz | grep "" -c` >> ${OUTPUT}/${j}.align.log

   ## Remove rRNA and reverse the strand (PRO-seq).
   gunzip -c ${TMPDIR}/$j.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${TMPDIR}/$j.nr.rs.bed.gz
   echo 'Number of mappable reads (excluding rRNA):' >> ${OUTPUT}/${j}.align.log
   echo `gunzip -c ${TMPDIR}/$j.nr.rs.bed.gz | grep "" -c` >> ${OUTPUT}/${j}.align.log

   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand + > ${TMPDIR}/$j\_plus.bedGraph
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand - > ${TMPDIR}/$j\_minus.noinv.bedGraph

   ## Invert minus strand.
   cat ${TMPDIR}/$j\_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${TMPDIR}/$j\_minus.bedGraph ## Invert read counts on the minus strand.

   ## Then to bigWig
   bedGraphToBigWig ${TMPDIR}/$j\_plus.bedGraph ${CHINFO} ${OUTPUT}/$j\_plus.bw
   bedGraphToBigWig ${TMPDIR}/$j\_minus.bedGraph ${CHINFO} ${OUTPUT}/$j\_minus.bw

  # rm ${TMPDIR}/$j.nr.rs.bed.gz ${TMPDIR}/$j.bed.gz ${TMPDIR}/$j*.bedGraph

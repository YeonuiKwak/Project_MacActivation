#!/usr/bin/bash
#
while test $# -gt 0; do
        case "$1" in
                -h|--help)
    echo ""
    echo "Preprocesses and aligns TED-seq data."
    echo "usage: TED-seq Preprocessig.sh -SE -genome ~/Work/shared/ref/hg38/ -c ~/Work/shared/ref/hg38/chrNameLength.txt -fastq ./*.fastq.gz -TED"
    echo "Takes PREFIX.fastq.gz (SE)"
    echo "or *.fastq.gz in the current working directory as input and writes"
    echo "BAM and bigWig files as output to the user_Assigned output directory"
    echo "requirement in the current working directory: cutadapt 1.8.3, fastx_trimmer, seqtk, prinseq-lite.pl 0.20.2, bwa, samtools, bedtools, and bedGraphToBigWig."
    echo "set PATH =$PATH:/home/labuser/anaconda2/bin/"
    echo "--UMI1=8 [The length of UMI barcode on the 5' end of R1 Read]"


                    exit 0
                    ;;
            -SE)
                    export SEQ="SE"
                    shift
                    ;;
            -PE)
                    export SEQ="PE"
                    shift
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
                            echo "no chromInfo specified; twp fieds;one for chromosome, the other for length info."
                            exit 1
                    fi
                    shift
                    ;;
            -fastq)
                    shift
                    if test $# -gt 0; then
                            export FQ_INPUT=$1
                    else
                            echo "fastq file is not found."
                            exit 1
                    fi
                    shift
                    ;;
            -TED)
                    export SE_OUTPUT="TED"
                    export RNA5="R1_5prime"
                    export OPP="FALSE"
                    export MAP5="TRUE"
                    export SE_READ="RNA_5prime"
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
    echo "--c is required"
    exit 1
fi
if [ -z "$SEQ" ]; then
    echo "please specify the input data is SE or PE"
    exit 1
fi
if [ -z "$SE_OUTPUT" ]; then
    echo "specify -TED option."
    exit 1
fi


#Preprocess TED-seq data
if [ -z "$UMI" ]; then
	export UMI=8
fi

echo " "
echo "processing PRO-seq data ..."
echo " "
echo "SEQ         $SEQ"
#SE
if [[ "$SEQ"=="SE" ]] ; then
    echo "SE_OUTPUT                 $SE_OUTPUT"
    echo "SE_READ                  $SE_READ"
fi
echo "Report 5' ends $MAP5"
echo "report opposite strand $OPP"
echo ""
echo "input/file paths:"
echo "star indexed genome directory         $STARIDX"
echo "chromInfo                            $CHINFO"

#Make directories
TMPDIR="./tmp"
mkdir ${TMPDIR}

mkdir ${TMPDIR}/noadapt
mkdir ${TMPDIR}/passQC

dedup_l=15
mkdir ${TMPDIR}/noadapt/l${dedup_l}
mkdir ${TMPDIR}/noadapt/l${dedup_l}_nodups
Name=`echo ${FQ_INPUT}| rev| cut -d\/ -f 1| rev`
mkdir ./${Name}_output

#########################################

echo ""
echo "preprocessing fastq files"
echo "trimmming poly(A) of read 1 > 10 nt"
##read stats
echo "${FQ_INPUT}" >${TMPDIR}/${Name}.QC.log
echo "Number of original input reads:" >>${TMPDIR}/${Name}.QC.log
gunzip -c ${FQ_INPUT}.fastq.gz| grep @ -c >>${TMPDIR}/${Name}.QC.log

##remove poly(A) and keep read length >=15 after trimming.
#There is no chance that the adaptor is read in TED-seq due to the size selection.#skipping adaptor removal
#poly(A) trimming #low-quality trimming??
gunzip -c ${FQ_INPUT}.fastq.gz >${FQ_INPUT}.fastq &
wait
prinseq-lite.pl -trim_tail_right 10 -fastq ${FQ_INPUT}.fastq -out_format 3 -out_bad first_bad -out_good ${TMPDIR}/noadapt/${Name}_trimmed_1 2> ./${Name}_output/prinseq-pcrdups.gd &
wait
echo "Number of reads after poly(A) trimming: " >>${TMPDIR}/${Name}.QC.log
cat ${TMPDIR}/noadapt/${Name}_trimmed_1.fastq|grep @ -c >>${TMPDIR}/${Name}.QC.log
echo ""
echo "filtering out reads with quality score >=20 and min length 20 nt"
prinseq-lite.pl -min_len 15 -min_qual_mean 20 -fastq ${TMPDIR}/noadapt/${Name}_trimmed_1.fastq -out_format 3 -out_bad ${TMPDIR}/noadapt/${Name}_trimmed_bad -out_good ${TMPDIR}/noadapt/${Name}_trimmed 2>> ./${Name}_output/prinseq-pcrdups.gd &
wait
# trimmed.fastq in the noadapt folder.
echo "Number of reads after poly(A) trimming and QC: " >>${TMPDIR}/${Name}.QC.log
cat ${TMPDIR}/noadapt/${Name}_trimmed.fastq|grep @ -c >>${TMPDIR}/${Name}.QC.log
########################################################################

#Collapse reads using prinseq-lite.pl #remove PCR duplicates.
echo "removing PCR duplciates started."
echo "extract first 15 nt."
cat ${TMPDIR}/noadapt/${Name}_trimmed.fastq |fastx_trimmer -Q33 -l ${dedup_l} -o ${TMPDIR}/noadapt/l${dedup_l}/${Name}_trimmed_l${dedup_l}.fastq &
wait
#deduplicate using the first dedup_l nt.
echo "removing PCR replicates using prinseq."
prinseq-lite.pl -fastq ${TMPDIR}/noadapt/l${dedup_l}/${Name}_trimmed_l${dedup_l}.fastq -derep 1 \
-out_format 3 -out_bad ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}_bad  -out_good ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l} 2>> ./${Name}_output/prinseq-pcrdups.gd &
wait
echo "deduplicate using the first 15 nt completed."
#dedup_l15.fastq
rm ${TMPDIR}/noadapt/l${dedup_l}/${Name}_trimmed_l${dedup_l}.fastq
#make a list of name from deduplicated fastq
echo "save sequence ids after PCR duplicates removal."
cat ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}.fastq \
|awk '(NR%4==1){print substr($1,2)}' > ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}.txt &
wait
#rm ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}.fastq
#generate fastq from the list of name
seqtk subseq ${TMPDIR}/noadapt/${Name}_trimmed.fastq ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}.txt > ${TMPDIR}/passQC/${Name}_dedup_withbarcode.fastq &
wait
echo 'Number of reads after PCR duplicates removal and QC:' >> ${TMPDIR}/${Name}.QC.log
cat ${TMPDIR}/passQC/${Name}_dedup_withbarcode.fastq | grep @ -c >> ${TMPDIR}/${Name}.QC.log

rm ${TMPDIR}/noadapt/l${dedup_l}_nodups/${Name}_dedup_l${dedup_l}.txt
rm ${TMPDIR}/noadapt/${Name}_trimmed.fastq

#./passQC/dedup_withbarcode.fastq
#trim the UMI from the 3' end of read1 after deduplicating.
echo "trimming UMI started."
prinseq-lite.pl -trim_left ${UMI} -fastq ${TMPDIR}/passQC/${Name}_dedup_withbarcode.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${Name}_dedup_barcoderemoved 2>> ./${Name}_output/prinseq-pcrdups.gd &
wait
#dedup_barcoderemoved.fastq
#minimum 15
prinseq-lite.pl -min_len 15 -fastq ${TMPDIR}/passQC/${Name}_dedup_barcoderemoved.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${Name}_dedup_QC_end 2>> ./${Name}_output/prinseq-pcrdups.gd &
wait
#dedup_QC_end.fastq
#rm ${TMPDIR}/passQC/${Name}_dedup_withbarcode.fastq
#rm ${TMPDIR}/passQC/${Name}_dedup_barcoderemoved.fastq

echo 'Number of reads after PCR duplicates, UMI removal and min 15 and QC:' >> ${TMPDIR}/${Name}.QC.log
cat ${TMPDIR}/passQC/${Name}_dedup_QC_end.fastq | grep @ -c >> ${TMPDIR}/${Name}.QC.log

#Cleanup
rm -r ${TMPDIR}/noadapt
#gzip ${TMPDIR}/passQC/${Name}_dedup_QC_end.fastq

###############Align reads.##########################
echo "Mapping reads:"
#align using star
~/Work/bin/star --genomeDir $STARIDX \
        --readFilesIn ${TMPDIR}/passQC/${Name}_dedup_QC_end.fastq \
        --outFilterMultimapNmax 1 --runThreadN 8 \
        --outFileNamePrefix ${TMPDIR}/passQC/${Name}. &
        wait
samtools view -b -q 10 ${TMPDIR}/passQC/${Name}.Aligned.out.sam|samtools sort -@ 8 - >${TMPDIR}/${Name}.sort.bam &
wait
########################
rm ${TMPDIR}/passQC/${Name}.Aligned.out.sam
cp ${TMPDIR}/${Name}.sort.bam ./${Name}_output/ &
wait
#cleanup
find ${TMPDIR} -name "*.sort.bam" -size -1024k -delete
#Bigwig#
#everything will be run in ${Name}_output
## Write out the bigWigs.
echo " "
echo "Writing bigWigs:"
bedtools bamtobed -i ${TMPDIR}/${Name}.sort.bam 2> ${TMPDIR}/kill.warnings|awk 'BEGIN{OFS="\t"}($5 > 0){print $0}'|awk 'BEGIN{OFS="\t"}($6 == "+"){print $1,$2,$2+1,$4,$5,$6;next}($6 == "-"){print $1,$3-1,$3,$4,$5,$6}'|gzip > ${TMPDIR}/${Name}_bed.gz
echo "number of reads :" >> ./${Name}_output/Align.log
echo `gunzip -c ${TMPDIR}/${Name}_bed.gz |grep "" -c` >>./${Name}_output/Align.log
## Remove rRNA and reverse the strand (PRO-seq).
gunzip -c ${TMPDIR}/${Name}_bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ./${Name}_output/${Name}_nr.rs.bed.gz
echo 'Number of mappable reads (excluding rRNA):'  >> ./${Name}_output/Align.log
echo `gunzip -c ./${Name}_output/${Name}_nr.rs.bed.gz | grep "" -c` >>./${Name}_output/Align.log
## Convert to bedGraph ... Can't gzip these, unfortunately.
bedtools genomecov -bg -i ./${Name}_output/${Name}_nr.rs.bed.gz  -g ~/Work/shared/ref/hg38/chrNameLength.txt -strand + > ./${Name}_output/${Name}.plus.bedGraph
bedtools genomecov -bg -i ./${Name}_output/${Name}_nr.rs.bed.gz  -g ~/Work/shared/ref/hg38/chrNameLength.txt -strand - > ./${Name}_output/${Name}.minus.noinv.bedGraph
## Invert minus strand.
cat ./${Name}_output/${Name}.minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,-1*$4}' > ./${Name}_output/${Name}.minus.bedGraph
## Invert readcounts on the minus strand.
## Then to bigWig.
#export CHINFO="~/Work/shared/ref/hg38/chrNameLength.txt"
#CHINFO="~/Work/shared/ref/hg38TEDi500/refMrna.500"
bedGraphToBigWig ./${Name}_output/${Name}.plus.bedGraph ${CHINFO} ./${Name}_output/${Name}_plus.bw
bedGraphToBigWig ./${Name}_output/${Name}.minus.bedGraph ${CHINFO} ./${Name}_output/${Name}_minus.bw
#clean up
rm ./${Name}_output/${Name}_nr.rs.bed.gz
#"rm ${TMPDIR}/${Name}_bed.gz"#

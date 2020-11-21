#!/usr/bin/bash
#
while test $# -gt 0; do
        case "$1" in
                -h|--help)
    echo ""
    echo "Calculte TED-seq mapped file(bam file) "
    echo "usage: TEDseq_PAL_estimation.sh -bam <##.bam> -ref <hg38.bed12>"
    echo "bash TEDseq_PAL_estimation.sh -bam <bam file> -ref ~/Work/labuser/sl2665/annotation/GencodeComprehensiveV26-hg38.bed"
    echo "sorted bam file will be used as uan input and hg38.bed12 will be used to create last 500 bp of transcripts based on 3'CPS annotation publicly available."
    echo ""
    echo "set PATH =$PATH:/home/labuser/anaconda2/bin/:~/Work/labuser/sl2665/stoat/bin/"
                    exit 0
                    ;;
            -bam)
                    shift
                    if test $# -gt 0; then
                            export INPUT=$1
                    else
                            echo "no input bam file is specified"
                            exit 1
                    fi
                    shift
                    ;;
            -ref)
                    shift
                    if test $# -gt 0; then
                            export REF=$1
                    else
                            echo "no hg38_annotation bed12 file specified"
                            exit 1
                    fi
                    shift
                    ;;
        esac
done
#Check arguments
if [ -z "$INPUT" ]; then
    echo "bam file is required."
    exit 1
fi
if [ -z "$REF" ]; then
    echo "ref annotation bed 12 is required."
    exit 1
fi


# Calculate PAL based on 3' CPS annotation-based last 500 bp window.
echo " "
echo "last window 500 based on 3'CPS annotation"

echo "input $INPUT"
echo "reference annotation $REF"
echo ""
#Make directories
TMPDIR="./tmp"
mkdir ${TMPDIR}

Name=`echo ${INPUT}| rev| cut -d\/ -f 1| rev`
mkdir ./${Name}_output


#######################################################\
#Sort bam file  by read name which means ordering by chromosome name.
samtools view ${INPUT} -H >tmp.bam
samtools view ${INPUT}|sort -k3,3 -k4,4n >>tmp.bam
bash ./tedseq-makepal3.sh -a tmp.bam -b ${REF} >./${Name}_output/tedseq.pal.txt

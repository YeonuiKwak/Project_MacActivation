#!/usr/bin/bash
#
while test $# -gt 0; do
        case "$1" in
                -h|--help)
    echo ""
    echo "Preprocesses and aligns PRO-seq data."
    echo ""
    echo "Takes PREFIX.fastq.gz (SE),  PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz (PE)"
    echo "or *.fastq.gz in the current working directory as input and writes" 
    echo "BAM and bigWig files as output to the user-assigned output-dir."
    echo ""
    echo "Requirements in current working directory:"
    echo "cutadapt 1.8.3, fastx_trimmer, seqtk, prinseq-lite.pl 0.20.2, bwa, samtools, bedtools, and bedGraphToBigWig."
    echo ""
    echo "bash proseq2.0.bsh [options]"
    echo ""
    echo "options:"
    echo ""
    echo "To get help:"
    echo "-h, --help             Show this brief help menu."
    echo ""
    echo "Required options:"
    echo "-SE, --SEQ=SE          Single-end sequencing."
    echo "-PE, --SEQ=PE          Paired-end sequencing."
    echo "-i, --bwa-index=PATH   Path to the BWA index of the target genome"
    echo "                       (i.e., bwa index)."
    echo "-c, --chrom-info=PATH  Location of the chromInfo table."
    echo ""
    echo "I/O options:"
    echo "-I, --fastq=PREFIX     Prefix for input files."
    echo "                       Paired-end files require identical prefix"
    echo "                       and end with _R1.fastq.gz and _R2.fastq.gz"
    echo "                       eg: PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz."
    echo "-T, --tmp=PATH         Path to a temporary storage directory."
    echo "-O, --output-dir=DIR   Specify a directory to store output in."
    echo ""
    echo "Required options for SE"
    echo "-G, --SE_READ=RNA_5prime Single-end sequencing from 5' end of"
    echo "                         nascent RNA, like GRO-seq."
    echo "-P, --SE_READ=RNA_3prime Single-end sequencing from 3' end of"
    echo "                         nascent RNA, like PRO-seq."
    echo ""
    echo "Options for PE"
    echo "--RNA5=R1_5prime    Specify the location of the 5' end of RNA"
    echo "                    [default: R1_5prime]."
    echo "--RNA3=R2_5prime    Specify the location of the 3' end of RNA"
    echo "                    [default: R2_5prime]."
    echo "                    Available options: R1_5prime: the 5' end of R1 reads"
    echo "                                       R2_5prime: the 5' end of R2 reads"
    echo "-5, --map5=TRUE     Report the 5' end of RNA [default on, --map5=TRUE]."
    echo "-3, --map5=FALSE    Report the 3' end of RNA,"
    echo "                    only available for PE [default off, --map5=TRUE]."
    echo "-s, --opposite-strand=TRUE"
    echo "                    Enable this option if the RNA are at the different strand"
    echo "                    as the reads set at RNA5 [default: disable]."
    echo ""
    echo "Optional operations:"
    echo "--ADAPT_SE=TGGAATTCTCGGGTGCCAAGG"
    echo "                    3' adapter to be removed from the 3' end of SE reads."
    echo "                   [default:TGGAATTCTCGGGTGCCAAGG]"
    echo "--ADAPT1=GATCGTCGGACTGTAGAACTCTGAACG"
    echo "                    3' adapter to be removed from the 3' end of R2."
    echo "                   [default:GATCGTCGGACTGTAGAACTCTGAACG]"
    echo "--ADAPT2=AGATCGGAAGAGCACACGTCTGAACTC"
    echo "                    3' adapter to be removed from the 3' end of R1."
    echo "                   [default:AGATCGGAAGAGCACACGTCTGAACTC]"
    echo ""
    echo "--UMI1=0           The length of UMI barcode on the 5' of R1 read. "
    echo "                   [default: 0]"
    echo "--UMI2=0           The length of UMI barcode on the 5' of R2 read. "
    echo "                   [default: 0]"
    echo "When UMI1 or UMI2 are set > 0, the pipeline will perform PCR deduplicate."
    echo ""
    echo "--Force_deduplicate=FALSE"
    echo "                   When --Force_deduplicate=TRUE, it will force the pipeline to"
    echo "                   perform PCR deduplicate even there is no UMI barcode"
    echo "                   (i.e. UMI1=0 and UMI2=0). [default: FALSE]"      
    echo "--ADD_B1=0         The length of additional barcode that will be trimmed"
    echo "                   on the 5' of R1 read. [default: 0]"
    echo "--ADD_B2=0         The length of additional barcode that will be trimmed"
    echo "                   on the 5' of R2 read. [default: 0]"
    echo "--thread=1         Number of threads can be used [default: 1]"
    echo ""
    echo "-4DREG             Using the pre-defined parameters to get the most reads"
    echo "                   for dREG package. Please use this flag to make the bigWig"
    echo "                   files compatible with dREG algorithm. [default: off]" 
    echo "-aln               Use BWA-backtrack [default: SE uses BWA-backtrack, PE uses BWA-MEM]"
    echo "-mem               Use BWA-MEM [default: SE uses BWA-backtrack, PE uses BWA-MEM]"



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
                -i)
                        shift
                        if test $# -gt 0; then
                                export BWAIDX=$1
                        else
                                echo "no BWA index specified"
                                exit 1
                        fi
                        shift
                        ;;
                --bwa-index*)
                        export BWAIDX=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -c)
                        shift
                        if test $# -gt 0; then
                                export CHINFO=$1
                        else
                                echo "no chromInfo specified"
                                exit 1
                        fi
                        shift
                        ;;
                --chrom-info*)
                        export CHINFO=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -I)
                        shift
                        if test $# -gt 0; then
                                export FQ_INPUT=$1
                        else
                                echo "no input prefix specified."
                                exit 1
                        fi
                        shift
                        ;;
                --fastq*)
                        export FQ_INPUT=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -T)
                        shift
                        if test $# -gt 0; then
                                export TMPDIR=$1
                        else
                                echo "no temp folder specified."
                                exit 1
                        fi
                        shift
                        ;;
                --tmp*)
                        export TMPDIR=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -O)
                        shift
                        if test $# -gt 0; then
                                export OUTPUT=$1
                        else
                                echo "no output dir specified."
                                exit 1
                        fi
                        shift
                        ;;
                --output-dir*)
                        export OUTPUT=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --thread*)
                        export thread=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --RNA5*) # report location of the 5 prime end of RNA
                         # acce
                        export RNA5=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --RNA3*) # report location of the 5 prime end of RNA
                         # acce
                        export RNA3=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -5)
                        export MAP5="TRUE"
                        shift
                        ;;
                --map5*)
                        export MAP5=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -3)
                        export MAP5="FALSE"
                        shift
                        ;;
                -s)
                        export OPP="TRUE"
                        shift
                        ;;
                --opposite-strand*)
                        export OPP=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --ADAPT_SE*)
                        export ADAPT_SE=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --ADAPT1*)
                        export ADAPT1=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --ADAPT2*)
                        export ADAPT2=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --UMI1*)
                        export UMI1=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --UMI2*)
                        export UMI2=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --Force_deduplicate*)
                        export Force_deduplicate=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                --ADD_B1*)
                        export ADD_B1=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;                
                --ADD_B2*)
                        export ADD_B2=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -G)
                        export SE_OUTPUT="G"
                        export RNA5="R1_5prime"
                        export OPP="FALSE"
                        #export MAP5="TRUE"
                        export SE_READ="RNA_5prime"
                        shift
                        ;;
                -P)
                        export SE_OUTPUT="P"
                        export RNA3="R1_5prime"
                        export OPP="TRUE"
                        #export MAP5="TRUE"
                        export SE_READ="RNA_3prime"
                        shift
                        ;;
                --SE_READ*)
                        export SE_READ=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -4DREG*)
                        export DREG_MODEL=1
                        shift
                        ;;
                -mem*)
                        export mem=1
                        shift
                        ;;
                -aln*)
                        export aln=1
                        shift
                        ;;
                *)
                        break
                        ;;
        esac
done


## CHECK ARGUMENTS.
if [ -z "$BWAIDX" ]; then
	echo "--bwa-index is required."
	echo " use bash proseq2.0.bsh --help."
	exit 1
fi
if [ -z "$CHINFO" ]; then
        echo "--chrom-info is required."
        echo " use bash proseq2.0.bsh --help."
        exit 1
fi
if [ -z "$SEQ" ]; then
  echo "Please specify the input data is SE or PE."
  echo " use bash proseq2.0.bsh --help."
  exit 1
fi

## INPUT & Parameters
# PE
if [[ "$SEQ" == "PE" ]] ; then 
  if [ -z "$SE_OUTPUT" ]; then 
    :
  else
    echo "-G and -P can only be used with -SE."
    echo " use bash proseq2.0.bsh --help."
    exit 1
  fi


  if [[ "$FQ_INPUT" == "*.fastq.gz" ]]; then
    FQ_INPUT=`ls *.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3-| cut -d _ -f 2- |rev | sort | uniq`
  fi
  if [ -z "$FQ_INPUT" ]; then
    echo "No input files specified.  Using *.fastq.gz"
    FQ_INPUT=`ls *.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3-| cut -d _ -f 2- |rev | sort | uniq`
  fi

  ## Check if input file end with _R1.fastq.gz or _R2.fastq.gz
  tmp=${FQ_INPUT}
  FQ_INPUT=() #array
  for name in ${tmp}
  do 
    if [ ! -f ${name}_R1.fastq.gz ]; then
      echo ""
      echo "##################################################################################"
      echo " File ${name}_R1.fastq.gz was not found! Skip ${name}*.fastq.gz from analysis."
      echo " Paired-end files require identical prefix and end with _R1.fastq.gz and _R2.fastq.gz"
      echo " eg: PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz"
      echo " Please make sure you have the correct file suffix!                           " 
      echo "##################################################################################"
      echo ""
    else
    FQ_INPUT+=(${name}) # append to the array
    #echo ${FQ_INPUT} 
    fi
  done

# SE
elif [[ "$SEQ" == "SE" ]] ; then
  if [[ "$SE_READ" == "RNA_5prime" ]] ; then
    SE_OUTPUT="G"
    RNA5="R1_5prime"
    OPP="FALSE"
    #MAP5="TRUE"
  elif [[ "$SE_READ" == "RNA_3prime" ]] ; then
    SE_OUTPUT="P"
    RNA3="R1_5prime"
    OPP="TRUE"
    #MAP5="TRUE"
  fi
  if [ -z "$SE_OUTPUT" ] ; then
  echo "Please specify output format for SE [-G or -P]"
  exit 1
  fi
  if [[ "$FQ_INPUT" == "*.fastq.gz" ]]; then
    FQ_INPUT=`ls *.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3-| rev`
  elif [ -z "$FQ_INPUT" ]; then
    echo "No input files specified.  Using *.fastq.gz"
    FQ_INPUT=`ls *.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3-| rev`
  elif [[ "$FQ_INPUT" == *.fastq.gz ]]; then
    FQ_INPUT=`echo $FQ_INPUT | rev | cut -d \. -f 3-| rev`
  fi
else
  echo "Please specify the input data is SE or PE."
  echo " use bash proseq2.0.bsh --help."
  exit 1
fi


# Check input file number
if [[ ${#FQ_INPUT[@]} == 0 ]]; then  # if the length of array is 0
  echo "##################################################################################"
  echo " No files is in the correct format."
  echo " Paired-end files require identical prefix and end with _R1.fastq.gz and _R2.fastq.gz"
  echo " eg: PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz"
  echo " Please make sure you have the correct file suffix! Aborting."
  echo "##################################################################################"
  exit 1
fi

if [ -z "$OUTPUT" ]; then
  now=$(date +"%m_%d_%Y")
  OUTPUT=./My_proseq_output_dir-${now}
  echo No output path specified.  Using ${OUTPUT}
fi

if [ ! -d $OUTPUT ]; then
  mkdir $OUTPUT
fi
if [ -z "$TMPDIR" ]; then
        TMPDIR="./"
fi
if [ ! -d $TMPDIR ]; then
  mkdir $TMPDIR
fi
# bash generate random 32 character alphanumeric string (upper and lowercase).
tmp=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1`
TMPDIR=$TMPDIR/$tmp

if [ -z "$thread" ]; then
   thread=1
fi

if [ -z "$ADAPT_SE" ]; then
  ADAPT_SE="TGGAATTCTCGGGTGCCAAGG"
fi
if [ -z "$ADAPT2" ]; then
  ADAPT2="AGATCGGAAGAGCACACGTCTGAACTC"
fi
if [ -z "$ADAPT1" ]; then
  ADAPT1="GATCGTCGGACTGTAGAACTCTGAACG"
fi
if [ -z "$UMI1" ]; then
   UMI1=0
elif [[ $[UMI1] -lt 0 ]]; then
  echo "UMI1 only take natural number"
  exit 1
fi
if [ -z "$UMI2" ]; then
   UMI2=0
elif [[ $[UMI2] -lt 0 ]]; then
  echo "UMI2 only take natural number"
  exit 1
fi
if [ -z "$Force_deduplicate" ]; then
   Force_deduplicate="FALSE"
elif [[ "$Force_deduplicate" == "TRUE" || "$Force_deduplicate" == "FALSE" ]] ; then
:
else
  echo "Force_deduplicate only take TRUE or FALSE"
  exit 1
fi


if [ -z "$ADD_B1" ]; then
   ADD_B1=0
fi
if [ -z "$ADD_B2" ]; then
   ADD_B2=0
fi


if [[ "$SEQ" == "PE" ]]; then
  if [[ -z "$RNA5" && -z "$RNA3" ]]; then
    RNA5="R1_5prime"
  fi
  if [[ "$RNA3" == "R1_5prime" ]]; then
   RNA5="R2_5prime"
  elif [[ "$RNA3" == "R2_5prime" ]]; then
    RNA5="R1_5prime"
  fi
  
  if [[ "$RNA5" == "R1_5prime" || "$RNA5" == "R2_5prime" ]] ; then 
          :
  else
          echo "--RNA5 and --RNA3 value can only be R1_5prime or R2_5prime."
          echo " use bash proseq2.0.bsh --help."
          exit 1
  fi
fi

if [ -z "$OPP" ]; then
   OPP="FALSE"
fi
if [ -z "$MAP5" ]; then
        MAP5="TRUE"
fi


if [[ "${MAP5}" == "FALSE" && "$SEQ" == "SE" ]] ; then 
  echo "For single-end (SE), can only report the 5 prime end of reads (--map5=TRUE)"
  exit 1
fi


if [ "${MAP5}" == "TRUE" ] ; then 
        :
elif [ "${MAP5}" == "FALSE" ] ; then
        :
else
        echo "--map5 value can only be TRUE or FALSE."
        echo " use bash proseq2.0.bsh --help."
        exit 1
fi

## Check all the bioinformatics tools can be called from current working directory.
for tool in cutadapt fastx_trimmer seqtk prinseq-lite bwa samtools bedtools bedGraphToBigWig sort-bed
do command -v ${tool} >/dev/null 2>&1 || { echo >&2 ${tool}" is required. Please make sure you can call the bioinformatics tools from your current working directoryb.  Aborting."; exit 1; }
done

#exec 1>test.log 2>&1
exec > >(tee ${OUTPUT}/proseq2.0_Run_${tmp}.log)
exec 2>&1
## Print out
echo " " 
echo "Processing PRO-seq data ..." 
echo " "

if [[ ! -z "$DREG_MODEL" ]]; then
    echo "dREG compatible mode      Yes"
fi

echo "SEQ                       $SEQ"

if [[ "$SEQ" == "SE" ]]; then
echo "SE_OUTPUT                 $SE_OUTPUT"
echo "SE_READ                   $SE_READ"
fi

if [[ "$SEQ" == "PE" ]]; then
echo "Location of 5' of RNA     $RNA5"
echo "Location of 3' of RNA     $RNA3"
fi

echo "Report 5' ends            $MAP5"
echo "Report opposite strand    $OPP"

echo ""
echo "Input files/ paths:"
echo "bwa index                 $BWAIDX"
echo "chromInfo                 $CHINFO"
if [[ "$SEQ" == "SE" ]]; then
i=1
for name in ${FQ_INPUT[@]}
  do 
echo "input file $i              ${name}.fastq.gz"
  let "i++"
done
fi
if [[ "$SEQ" == "PE" ]]; then
i=1
for name in ${FQ_INPUT[@]}
  do 
echo "input file pair $i         ${name}_R1.fastq.gz, ${name}_R2.fastq.gz"
  let "i++"
done
fi
echo "temp folder               $TMPDIR"
echo "output-dir                $OUTPUT"
echo " "
echo "Optional operations:"
echo "ADAPT_SE                  $ADAPT_SE"
echo "ADAPT1                    $ADAPT1"
echo "ADAPT2                    $ADAPT2"
echo "UMI1 barcode length       $UMI1"
echo "UMI2 barcode length       $UMI2"
echo "ADD_B1 length             $ADD_B1"
echo "ADD_B2 length             $ADD_B2"
echo "number of threads         $thread"
if [[ ${UMI2} != 0 || ${UMI1} != 0 || ${Force_deduplicate} == "TRUE" ]]; then 
 echo "Remove PCR duplicates     TRUE"
else
 echo "Remove PCR duplicates     FALSE"
fi


## Exits .. for debugging.
#exit 1

## DOIT!
mkdir ${TMPDIR}

if [[ "$SEQ" == "PE" ]] ; then 
#############################################
## Preprocess data.  Remove adapters.  Trim.
echo " "
echo "Preprocessing fastq files:"
mkdir ${TMPDIR}/noadapt
mkdir ${TMPDIR}/passQC
if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then 
  dedup_L=30 
  mkdir ${TMPDIR}/noadapt/l${dedup_L}
  mkdir ${TMPDIR}/noadapt/l${dedup_L}_nodups
fi

for name in ${FQ_INPUT[@]}
 do
  old_name=${name}
  name=`echo ${old_name} | rev | cut -d \/ -f 1 |rev`
  ## read stats
  echo 'Number of original input reads:' > ${TMPDIR}/${name}.QC.log
  echo ${old_name}_R1.fastq.gz >> ${TMPDIR}/${name}.QC.log
  zcat ${old_name}_R1.fastq.gz | grep @ -c >> ${TMPDIR}/${name}.QC.log
  echo ${old_name}_R2.fastq.gz >> ${TMPDIR}/${name}.QC.log
  zcat ${old_name}_R2.fastq.gz | grep @ -c >> ${TMPDIR}/${name}.QC.log



  ## Remove adapter, UMI barcode, additional barcode, and low quality (q=20) base from 3prime end of reads. Keep read length >=15 after trimmming
  # Remove adapter
  cutadapt -a ${ADAPT2} -e 0.10 --overlap 2 --output=${TMPDIR}/${name}_trim_R1.fastq --untrimmed-output=${TMPDIR}/${name}_untrim_R1.fastq ${old_name}_R1.fastq.gz &
  cutadapt -a ${ADAPT1} -e 0.10 --overlap 2 --output=${TMPDIR}/${name}_trim_R2.fastq --untrimmed-output=${TMPDIR}/${name}_untrim_R2.fastq ${old_name}_R2.fastq.gz &
  wait

  # Read1
  # remove UMI2 and ADD_B2 from the 3 prime end of R1
   n2=$[UMI2+ADD_B2]
   cutadapt --cut -${n2} --minimum-length=10 ${TMPDIR}/${name}_trim_R1.fastq --output=${TMPDIR}/${name}_trim.${n2}Nremoved_R1.fastq -q 20 &
   cutadapt --minimum-length=10 ${TMPDIR}/${name}_untrim_R1.fastq --output=${TMPDIR}/${name}_q20trim_R1.fastq -q 20 &
   wait
   cat ${TMPDIR}/${name}_q20trim_R1.fastq ${TMPDIR}/${name}_trim.${n2}Nremoved_R1.fastq | paste - - - - |LC_ALL=C sort --temporary-directory=${TMPDIR} --parallel=10 -k1,1 -S 10G | tr '\t' '\n' > ${TMPDIR}/noadapt/${name}_noadapt_R1.fastq &

  # Read2
  # remove UMI1 and ADD_B1 from the 3 prime end of R2
   n1=$[UMI1+ADD_B1]
   cutadapt --cut -${n1} --minimum-length=10 ${TMPDIR}/${name}_trim_R2.fastq --output=${TMPDIR}/${name}_trim.${n1}Nremoved_R2.fastq -q 20 &
   cutadapt --minimum-length=10 ${TMPDIR}/${name}_untrim_R2.fastq --output=${TMPDIR}/${name}_q20trim_R2.fastq -q 20 &
   wait
   cat ${TMPDIR}/${name}_q20trim_R2.fastq ${TMPDIR}/${name}_trim.${n1}Nremoved_R2.fastq | paste - - - - | LC_ALL=C sort --temporary-directory=${TMPDIR} --parallel=10 -k1,1 -S 10G | tr '\t' '\n' > ${TMPDIR}/noadapt/${name}_noadapt_R2.fastq &

  wait
  echo 'Number of reads after adapter removal and QC:' >> ${TMPDIR}/${name}.QC.log
  echo R1 >> ${TMPDIR}/${name}.QC.log
  cat ${TMPDIR}/noadapt/${name}_noadapt_R1.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log
  echo R2 >> ${TMPDIR}/${name}.QC.log
  cat ${TMPDIR}/noadapt/${name}_noadapt_R2.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log

  rm ${TMPDIR}/${name}_trim_R1.fastq ${TMPDIR}/${name}_untrim_R1.fastq ${TMPDIR}/${name}_q20trim_R1.fastq ${TMPDIR}/${name}_trim.${n2}Nremoved_R1.fastq &
  rm ${TMPDIR}/${name}_trim_R2.fastq ${TMPDIR}/${name}_untrim_R2.fastq ${TMPDIR}/${name}_q20trim_R2.fastq ${TMPDIR}/${name}_trim.${n1}Nremoved_R2.fastq &
  

  ## Collapse reads using prinseq-lite.pl. if there are UMI barcodes or $Force_deduplicate=TRUE
  if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then # if there is UMI barcode, Will perform deduplicates with first ${dedup_L} bp
   # Remove PCR duplciates.
   
   # fastq with the first dedup_L nt
   cat ${TMPDIR}/noadapt/${name}_noadapt_R1.fastq | fastx_trimmer -Q33 -l ${dedup_L} -o ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R1.fastq &
   cat ${TMPDIR}/noadapt/${name}_noadapt_R2.fastq | fastx_trimmer -Q33 -l ${dedup_L} -o ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R2.fastq &
   wait
   
   # deduplicate using the first dedup_L nt
   prinseq-lite -fastq ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R1.fastq -fastq2 ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R2.fastq -derep 1 -out_format 3 -out_bad null -out_good ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup -min_len 15 2> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   rm ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R1.fastq ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R2.fastq &
   # make a list of name from deduplicated fastq
   cat ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_1.fastq | awk '(NR%4==1){print substr($1,2)}' > ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt &
   cat ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_1_singletons.fastq | awk '(NR%4==1){print substr($1,2)}' > ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_1_singletons_l${dedup_L}.txt &
   cat ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_2_singletons.fastq | awk '(NR%4==1){print substr($1,2)}' > ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_2_singletons_l${dedup_L}.txt &
   wait
   rm ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_1.fastq ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_1_singletons.fastq ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_2_singletons.fastq

   # generate fastq from the list of name
    seqtk subseq ${TMPDIR}/noadapt/${name}_noadapt_R1.fastq ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt > ${TMPDIR}/passQC/${name}_dedup_withBarcode_1.fastq &
    seqtk subseq ${TMPDIR}/noadapt/${name}_noadapt_R2.fastq ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt > ${TMPDIR}/passQC/${name}_dedup_withBarcode_2.fastq &
    seqtk subseq ${TMPDIR}/noadapt/${name}_noadapt_R1.fastq ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_1_singletons_l${dedup_L}.txt > ${TMPDIR}/passQC/${name}_dedup_1_singletons.fastq &
    seqtk subseq ${TMPDIR}/noadapt/${name}_noadapt_R2.fastq ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_2_singletons_l${dedup_L}.txt > ${TMPDIR}/passQC/${name}_dedup_2_singletons.fastq &
    wait
    rm ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup*l${dedup_L}.txt &
     
    # trim the UMI and additional barcode after dereplicate
     prinseq-lite -trim_left ${n1} -fastq ${TMPDIR}/passQC/${name}_dedup_withBarcode_1.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_1 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd &
     prinseq-lite -trim_left ${n2} -fastq ${TMPDIR}/passQC/${name}_dedup_withBarcode_2.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_2 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd &
     wait
     rm ${TMPDIR}/passQC/${name}_dedup_withBarcode_1.fastq ${TMPDIR}/passQC/${name}_dedup_withBarcode_2.fastq
    # min_len 15
     prinseq-lite -min_len 15 -fastq ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_1.fastq -fastq2 ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_2.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_QC_end 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
     rm ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_1.fastq ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_2.fastq &
     echo 'Number of paired reads after PCR duplicates removal and QC:' >> ${TMPDIR}/${name}.QC.log
     cat ${TMPDIR}/passQC/${name}_dedup_QC_end_1.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log &

  elif [[ ${Force_deduplicate} == "TRUE" ]]; then # if there is NO UMI barcode, Will perform deduplicates with whole legnth of reads
   # Remove PCR duplciates.
   prinseq-lite -derep 1 -fastq ${TMPDIR}/noadapt/${name}_noadapt_R1.fastq -fastq2 ${TMPDIR}/noadapt/${name}_noadapt_R2.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_withBarcode 2> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   # trim the UMI and additional barcode after dereplicate
   prinseq-lite -trim_left ${n1} -fastq ${TMPDIR}/passQC/${name}_dedup_withBarcode_1.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_1 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   prinseq-lite -trim_left ${n2} -fastq ${TMPDIR}/passQC/${name}_dedup_withBarcode_2.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_2 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   rm ${TMPDIR}/passQC/${name}_dedup_withBarcode_1.fastq ${TMPDIR}/passQC/${name}_dedup_withBarcode_2.fastq
   # min_len 15
   prinseq-lite -min_len 15 -fastq ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_1.fastq -fastq2 ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_2.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_QC_end 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   rm ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_1.fastq ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved_2.fastq
   echo 'Number of paired reads after PCR duplicates removal and QC:' >> ${TMPDIR}/${name}.QC.log
   cat ${TMPDIR}/passQC/${name}_dedup_QC_end_1.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log
  else
   # trim the additional barcode ${ADD_B1} and ${ADD_B2} WITHOUT dereplicate. If no barcode, prinseq-lite.pl will remove unpair reads and reads that are length 0
   prinseq-lite -trim_left ${ADD_B1} -fastq ${TMPDIR}/noadapt/${name}_noadapt_R1.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_BarcodeRemoved_1 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   prinseq-lite -trim_left ${ADD_B2} -fastq ${TMPDIR}/noadapt/${name}_noadapt_R2.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_BarcodeRemoved_2 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   # min_len 15
   prinseq-lite -min_len 15 -fastq ${TMPDIR}/passQC/${name}_BarcodeRemoved_1.fastq -fastq2 ${TMPDIR}/passQC/${name}_BarcodeRemoved_2.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_QC_end 2> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   rm ${TMPDIR}/passQC/${name}_BarcodeRemoved_1.fastq ${TMPDIR}/passQC/${name}_BarcodeRemoved_2.fastq
   echo 'Number of paired reads after final QC:' >> ${TMPDIR}/${name}.QC.log
   cat ${TMPDIR}/passQC/${name}_QC_end_1.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log
  fi

  cat ${OUTPUT}/${name}.prinseq-pcrDups.gd

 done

wait
## Cleanup.
rm -r ${TMPDIR}/noadapt 
gzip ${TMPDIR}/passQC/*.fastq 


#############################################
## Align reads.
echo " "
echo "Mapping reads:"

QC_INPUT=`ls  ${TMPDIR}/passQC/*_QC_end_1.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- | cut -d _ -f 2- | rev| sort | uniq`

if [[ ${#QC_INPUT} == 0 ]]; then
  echo "#########################################"
  echo " Something went wrong with Preprocess."
  echo " No fastq.gz files passed."
  echo " Aborting."
  echo "#########################################"
  exit 1
fi

if [[ ! -z "$mem" ]]; then #if -mem was given, use mem
  for name in ${QC_INPUT}
    do
    ## Align using BWA.
    bwa mem -k 19 -t ${thread} ${BWAIDX} ${TMPDIR}/passQC/${name}_1.fastq.gz ${TMPDIR}/passQC/${name}_2.fastq.gz | \
    samtools view -bf 0x2 -q 20 - | samtools sort -n -@ ${thread} - > ${TMPDIR}/${name}.sort.bam
  done

elif [[ ! -z "$aln" ]]; then # elif -aln was given, use aln
 for name in ${QC_INPUT}
  do
  ## Align using BWA aln
  bwa aln -t ${thread} ${BWAIDX} ${TMPDIR}/passQC/${name}_1.fastq.gz  >  ${TMPDIR}/${name}_aln_sa1.sai
  bwa aln -t ${thread} ${BWAIDX} ${TMPDIR}/passQC/${name}_2.fastq.gz  >  ${TMPDIR}/${name}_aln_sa2.sai
  bwa sampe -n 1 -f ${TMPDIR}/passQC/${name}_end.sam ${BWAIDX} ${TMPDIR}/${name}_aln_sa1.sai ${TMPDIR}/${name}_aln_sa2.sai ${TMPDIR}/passQC/${name}_1.fastq.gz ${TMPDIR}/passQC/${name}_2.fastq.gz

  samtools view -b -q 20 ${TMPDIR}/passQC/${name}_end.sam | samtools sort -n -@ ${thread} - > ${TMPDIR}/${name}.sort.bam
  rm ${TMPDIR}/${name}_aln_sa1.sai ${TMPDIR}/${name}_aln_sa2.sai ${TMPDIR}/passQC/${name}_end.sam
  done

else #defaul use mem
  for name in ${QC_INPUT}
    do
    ## Align using BWA.
    bwa mem -k 19 -t ${thread} ${BWAIDX} ${TMPDIR}/passQC/${name}_1.fastq.gz ${TMPDIR}/passQC/${name}_2.fastq.gz | \
    samtools view -bf 0x2 -q 20 - | samtools sort -n -@ ${thread} - > ${TMPDIR}/${name}.sort.bam
  done
fi


for name in ${QC_INPUT}
  do
  cp ${TMPDIR}/${name}.sort.bam ${OUTPUT} &
done
wait







## Cleanup
find ${TMPDIR} -name "*.sort.bam" -size -1024k -delete

#############################################
## Write out the bigWigs.
echo " "
echo "Writing bigWigs:"
for f in ${TMPDIR}/*.sort.bam
 do
   j=`echo $f | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev `
   echo $j > ${OUTPUT}/${j}.align.log
   if [ "${RNA5}" == "R1_5prime" ] ; then
     if [ "${OPP}" == "FALSE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then  ## report The 5' end of the RNA. Danko lab leChRO-Seq protocol is on the 5' of _R1 readl, same strand of R1 ($9) 
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$9}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$9}' | gzip > ${TMPDIR}/$j.bed.gz
       else ## report The 3' end of the RNA.  Danko lab leChRO-Seq protocol is on the 5 prime of _R2 read, opposite strand of R2 (R2 strand $10, R1 strand $9)
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "-") {print $1,$6-1,$6,$7,$8,$9}; ($10 == "+") {print $1,$5,$5+1,$7,$8,$9}' | gzip > ${TMPDIR}/$j.bed.gz
       fi
     elif [ "${OPP}" == "TRUE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then  ## report The 5' end of the RNA. 
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$10}' | gzip > ${TMPDIR}/$j.bed.gz
       else ## report The 3' end of the RNA.
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "-") {print $1,$6-1,$6,$7,$8,$10}; ($10 == "+") {print $1,$5,$5+1,$7,$8,$10}' | gzip > ${TMPDIR}/$j.bed.gz
       fi
     fi
   elif [ "${RNA5}" == "R2_5prime" ] ; then
     if [ "${OPP}" == "FALSE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then #report the 5 prime end of RNA, in Engreitz data is 5 prime end of R2, same strand
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/${j}.kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$5,$5+1,$7,$8,$10}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$10}'|gzip > ${TMPDIR}/${j}.bed.gz      
       else  ## report the 3-prime end of the RNA, in Engreitz data is the 5' end of R1 read, but opposite strand
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/${j}.kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$10}' |gzip  > ${TMPDIR}/${j}.bed.gz
       fi
     elif [ "${OPP}" == "TRUE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then #report the 5 prime end of RNA, in Engreitz data is 5 prime end of R2, same strand
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/${j}.kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$5,$5+1,$7,$8,$9}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$9}'|gzip > ${TMPDIR}/${j}.bed.gz      
       else  ## report the 3-prime end of the RNA, in Engreitz data is the 5' end of R1 read, but opposite strand
          bedtools bamtobed -bedpe -mate1 -i $f 2> ${TMPDIR}/${j}.kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$9}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$9}' |gzip  > ${TMPDIR}/${j}.bed.gz
       fi
     fi
   fi

   echo 'Number of mappable reads:' >> ${OUTPUT}/${j}.align.log
   echo `zcat ${TMPDIR}/$j.bed.gz | grep "" -c` >> ${OUTPUT}/${j}.align.log
   
   ## Remove rRNA and reverse the strand (PRO-seq).
   zcat ${TMPDIR}/$j.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${TMPDIR}/$j.nr.rs.bed.gz
   echo 'Number of mappable reads (excluding rRNA):' >> ${OUTPUT}/${j}.align.log
   echo `zcat ${TMPDIR}/$j.nr.rs.bed.gz | grep "" -c` >> ${OUTPUT}/${j}.align.log
   
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand + > ${TMPDIR}/$j\_plus.bedGraph
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand - > ${TMPDIR}/$j\_minus.noinv.bedGraph
   
   ## Invert minus strand.
   cat ${TMPDIR}/$j\_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${TMPDIR}/$j\_minus.bedGraph ## Invert read counts on the minus strand.
   
   ## Then to bigWig
   bedGraphToBigWig ${TMPDIR}/$j\_plus.bedGraph ${CHINFO} ${OUTPUT}/$j\_plus.bw
   bedGraphToBigWig ${TMPDIR}/$j\_minus.bedGraph ${CHINFO} ${OUTPUT}/$j\_minus.bw

  #rm ${TMPDIR}/$j.nr.rs.bed.gz ${TMPDIR}/$j.bed.gz ${TMPDIR}/$j*.bedGraph
 done

fi #*/


if [[ "$SEQ" == "SE" ]] ; then 
#############################################
## Preprocess data.  Remove adapters.  Trim.
echo " "
echo "Preprocessing fastq files:"
mkdir ${TMPDIR}/noadapt
mkdir ${TMPDIR}/passQC
if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then 
  dedup_L=30 
  mkdir ${TMPDIR}/noadapt/l${dedup_L}
  mkdir ${TMPDIR}/noadapt/l${dedup_L}_nodups
fi

for name in ${FQ_INPUT[@]}
 do
  old_name=${name}
  name=`echo ${old_name} | rev | cut -d \/ -f 1| rev`
  ## read stats
  echo ${old_name} >  ${TMPDIR}/${name}.QC.log
  echo 'Number of original input reads:' >> ${TMPDIR}/${name}.QC.log
  zcat ${old_name}.fastq.gz | grep @ -c >> ${TMPDIR}/${name}.QC.log

  ## Remove adapter, UMI barcode, additional barcode, and low quality (q=20) base from 3prime end of reads. Keep read length >=15 after trimmming
  # Remove adapter
  cutadapt -a ${ADAPT_SE} -e 0.10 --overlap 2 --output=${TMPDIR}/${name}_trim.fastq --untrimmed-output=${TMPDIR}/${name}_untrim.fastq ${old_name}.fastq.gz
  # Read1
  # remove UMI2 and ADD_B2 from the 3 prime end of R1
   n2=$[UMI2+ADD_B2]
   n1=$[UMI1+ADD_B1]
   cutadapt --cut -${n2} --minimum-length=10 ${TMPDIR}/${name}_trim.fastq --output=${TMPDIR}/${name}_trim.${n2}Nremoved.fastq --nextseq-trim 20 &
   cutadapt --minimum-length=10 ${TMPDIR}/${name}_untrim.fastq --output=${TMPDIR}/${name}_q20trim.fastq --nextseq-trim 20 &
   wait
   cat ${TMPDIR}/${name}_q20trim.fastq ${TMPDIR}/${name}_trim.${n2}Nremoved.fastq | paste - - - - |LC_ALL=C sort --temporary-directory=${TMPDIR} --parallel=10 -k1,1 -S 10G | tr '\t' '\n' > ${TMPDIR}/noadapt/${name}_noadapt.fastq &

  wait
  echo 'Number of reads after adapter removal and QC:' >> ${TMPDIR}/${name}.QC.log
  cat ${TMPDIR}/noadapt/${name}_noadapt.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log


  rm ${TMPDIR}/${name}_trim.fastq ${TMPDIR}/${name}_untrim.fastq ${TMPDIR}/${name}_q20trim.fastq ${TMPDIR}/${name}_trim.${n2}Nremoved.fastq
  

  ## Collapse reads using prinseq-lite.pl. if there are UMI barcodes
  if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then 
   # Remove PCR duplciates.
   # fastq with the first dedup_L nt
   cat ${TMPDIR}/noadapt/${name}_noadapt.fastq | fastx_trimmer -Q33 -l ${dedup_L} -o ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}.fastq 
   # deduplicate using the first dedup_L nt
   prinseq-lite -fastq ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}.fastq -derep 1 -out_format 3 -out_bad null -out_good ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup -min_len 15 2> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   rm ${TMPDIR}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}.fastq

   # make a list of name from deduplicated fastq
   cat ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup.fastq | awk '(NR%4==1){print substr($1,2)}' > ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt
   rm ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup.fastq

   # generate fastq from the list of name
   seqtk subseq ${TMPDIR}/noadapt/${name}_noadapt.fastq ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt > ${TMPDIR}/passQC/${name}_dedup_withBarcode.fastq 
   rm ${TMPDIR}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt

   # trim the UMI and additional barcode after dereplicate
   prinseq-lite -trim_left ${n1} -fastq ${TMPDIR}/passQC/${name}_dedup_withBarcode.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   # min_len 15
   prinseq-lite -min_len 15 -fastq ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_QC_end 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   rm ${TMPDIR}/passQC/${name}_dedup_withBarcode.fastq ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved.fastq
   echo 'Number of reads after PCR duplicates removal and QC:' >> ${TMPDIR}/${name}.QC.log
   cat ${TMPDIR}/passQC/${name}_dedup_QC_end.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log
  
  elif [[ ${Force_deduplicate} == "TRUE" ]]; then 
   # Remove PCR duplciates.
   prinseq-lite -derep 1 -fastq ${TMPDIR}/noadapt/${name}_noadapt.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_withBarcode 2> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   # trim the UMI and additional barcode after dereplicate
   prinseq-lite -trim_left ${n1} -fastq ${TMPDIR}/passQC/${name}_dedup_withBarcode.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   # min_len 15
   prinseq-lite -min_len 15 -fastq ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_dedup_QC_end 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   rm ${TMPDIR}/passQC/${name}_dedup_withBarcode.fastq ${TMPDIR}/passQC/${name}_dedup_BarcodeRemoved.fastq
   echo 'Number of reads after PCR duplicates removal and QC:' >> ${TMPDIR}/${name}.QC.log
   cat ${TMPDIR}/passQC/${name}_dedup_QC_end.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log
  else
   # trim the additional barcode WITHOUT dereplicate. If no barcode, prinseq-lite.pl will remove unpair reads and reads that are length 0
   prinseq-lite -trim_left ${ADD_B1} -fastq ${TMPDIR}/noadapt/${name}_noadapt.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_BarcodeRemoved 2>> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   # min_len 15
   prinseq-lite -min_len 15 -fastq ${TMPDIR}/passQC/${name}_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${TMPDIR}/passQC/${name}_QC_end 2> ${OUTPUT}/${name}.prinseq-pcrDups.gd
   rm ${TMPDIR}/passQC/${name}_BarcodeRemoved.fastq
   echo 'Number of reads after final QC:' >> ${TMPDIR}/${name}.QC.log
   cat ${TMPDIR}/passQC/${name}_QC_end.fastq | grep @ -c >> ${TMPDIR}/${name}.QC.log
  fi

  cat ${OUTPUT}/${name}.prinseq-pcrDups.gd



 done

wait
## Cleanup.
rm -r ${TMPDIR}/noadapt 
gzip ${TMPDIR}/passQC/*.fastq 


#############################################
## Align reads.
echo " "
echo "Mapping reads:"

QC_INPUT=`ls  ${TMPDIR}/passQC/*_QC_end.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- | cut -d _ -f 2- | rev| sort | uniq`

if [[ ${#QC_INPUT} == 0 ]]; then
  echo "#########################################"
  echo " Something went wrong with Preprocess."
  echo " No fastq.gz files passed."
  echo " Aborting."
  echo "#########################################"
  exit 1
fi

if [[ -z "$DREG_MODEL" ]]; then
  if [[ ! -z "$aln" ]]; then #use aln
    for name in ${QC_INPUT}
      do
      ## Align using BWA aln
      bwa aln -t ${thread} ${BWAIDX} ${TMPDIR}/passQC/${name}_end.fastq.gz | \
      bwa samse -n 1 -f ${TMPDIR}/passQC/${name}_end.sam ${BWAIDX} - ${TMPDIR}/passQC/${name}_end.fastq.gz &
      done
    wait
    for name in ${QC_INPUT}
      do
      samtools view -b -q 20 ${TMPDIR}/passQC/${name}_end.sam |\
      samtools sort -n -@ ${thread} - > ${TMPDIR}/${name}.sort.bam &
      done
    wait
    for name in ${QC_INPUT}
      do
      rm ${TMPDIR}/passQC/${name}_end.sam
      cp ${TMPDIR}/${name}.sort.bam ${OUTPUT} & ## Saves the sorted BAM in the output file.
    done
    wait

  elif [[ ! -z "$mem" ]]; then #use mem
    for name in ${QC_INPUT}
      do
      ## Align using BWA.
      bwa mem -k 19 -t ${thread} ${BWAIDX} ${TMPDIR}/passQC/${name}_end.fastq.gz | \
      samtools view -b -q 20 - | samtools sort -n -@ ${thread} - > ${TMPDIR}/${name}.sort.bam &
      done
      wait
    for name in ${QC_INPUT}
      do
      cp ${TMPDIR}/${name}.sort.bam ${OUTPUT} & ## Saves the sorted BAM in the output file.
      done
      wait

  else #default use aln
    for name in ${QC_INPUT}
      do
      ## Align using BWA aln
      bwa aln -t ${thread} ${BWAIDX} ${TMPDIR}/passQC/${name}_end.fastq.gz | \
      bwa samse -n 1 -f ${TMPDIR}/passQC/${name}_end.sam ${BWAIDX} - ${TMPDIR}/passQC/${name}_end.fastq.gz &
      done
    wait
    for name in ${QC_INPUT}
      do
      samtools view -b -q 20 ${TMPDIR}/passQC/${name}_end.sam | samtools sort -n -@ ${thread} - > ${TMPDIR}/${name}.sort.bam &
      done
    wait
    for name in ${QC_INPUT}
      do
      rm ${TMPDIR}/passQC/${name}_end.sam
      cp ${TMPDIR}/${name}.sort.bam ${OUTPUT} & ## Saves the sorted BAM in the output file.
    done
    wait
 fi
else #dREG
  echo "Aligning reads using the same parameters with dREG" 
  for name in ${QC_INPUT}
    do
    ## Align using BWA aln , same as the original dreg model.
    bwa aln -t ${thread} ${BWAIDX} ${TMPDIR}/passQC/${name}_end.fastq.gz | \
    bwa samse -n 1 -f ${TMPDIR}/passQC/${name}_end.sam ${BWAIDX} - ${TMPDIR}/passQC/${name}_end.fastq.gz &
    done
  
  wait

  for name in ${QC_INPUT}
    do
    samtools view -b -q 0 ${TMPDIR}/passQC/${name}_end.sam | samtools sort -n -@ ${thread} - > ${TMPDIR}/${name}.sort.bam &
    done
  
  wait
  for name in ${QC_INPUT}
   do
   rm ${TMPDIR}/passQC/${name}_end.sam
   cp ${TMPDIR}/${name}.sort.bam ${OUTPUT} & ## Saves the sorted BAM in the output file.
   done
   wait
fi

wait


## Cleanup
find ${TMPDIR} -name "*.sort.bam" -size -1024k -delete


#############################################
## Write out the bigWigs.
echo " "
echo "Writing bigWigs:"
for f in ${TMPDIR}/*.sort.bam
 do
#*/
   j=`echo $f | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev `
   echo $j > ${OUTPUT}/${j}.align.log

# in SE, MAP5 alwasys TRUE


   #if [[ "${RNA5}" == "R1_5prime" && "${OPP}" == "FALSE" ]] ; then ## report The 5 prime end of the RNA.   #like GRO-seq
   if [[ "$SE_OUTPUT" == "G" ]] ; then
      bedtools bamtobed -i $f 2> ${TMPDIR}/kill.warnings| awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
      awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | gzip > ${TMPDIR}/$j.bed.gz
   #elif [[ "${RNA3}" == "R1_5prime" && "${OPP}" == "TRUE" ]] ; then  #like PRO-seq
    elif [[ "$SE_OUTPUT" == "P" ]] ; then
      bedtools bamtobed -i $f 2> ${TMPDIR}/kill.warnings| awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
      awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,"-"}; ($6 == "-") {print $1,$3-1,$3,$4,$5,"+"}' | gzip > ${TMPDIR}/$j.bed.gz
   fi


   echo 'Number of mappable reads:' >> ${OUTPUT}/${j}.align.log
   echo `zcat ${TMPDIR}/$j.bed.gz | grep "" -c` >> ${OUTPUT}/${j}.align.log
   
   ## Remove rRNA and reverse the strand (PRO-seq).
   zcat ${TMPDIR}/$j.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${TMPDIR}/$j.nr.rs.bed.gz
   echo 'Number of mappable reads (excluding rRNA):' >> ${OUTPUT}/${j}.align.log
   echo `zcat ${TMPDIR}/$j.nr.rs.bed.gz | grep "" -c` >> ${OUTPUT}/${j}.align.log
   
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand + > ${TMPDIR}/$j\_plus.bedGraph
   bedtools genomecov -bg -i ${TMPDIR}/$j.nr.rs.bed.gz -g ${CHINFO} -strand - > ${TMPDIR}/$j\_minus.noinv.bedGraph
   
   ## Invert minus strand.
   cat ${TMPDIR}/$j\_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${TMPDIR}/$j\_minus.bedGraph ## Invert read counts on the minus strand.
   
   ## Then to bigWig
   bedGraphToBigWig ${TMPDIR}/$j\_plus.bedGraph ${CHINFO} ${OUTPUT}/$j\_plus.bw
   bedGraphToBigWig ${TMPDIR}/$j\_minus.bedGraph ${CHINFO} ${OUTPUT}/$j\_minus.bw

  # rm ${TMPDIR}/$j.nr.rs.bed.gz ${TMPDIR}/$j.bed.gz ${TMPDIR}/$j*.bedGraph
 done



fi

echo "QC" > ${OUTPUT}/proseq2.0_read_report_${tmp}.log
for old_name in ${FQ_INPUT[@]}
 do
  name=`echo ${old_name} | rev | cut -d \/ -f 1 |rev`
  cat ${TMPDIR}/${name}.QC.log >> ${OUTPUT}/proseq2.0_read_report_${tmp}.log
 done
echo "" >> ${OUTPUT}/proseq2.0_read_report_${tmp}.log
echo "Mapping" >> ${OUTPUT}/proseq2.0_read_report_${tmp}.log
for f in ${TMPDIR}/*.sort.bam
 do j=`echo $f | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev `
 cat ${OUTPUT}/${j}.align.log >> ${OUTPUT}/proseq2.0_read_report_${tmp}.log
 done



#############################################
## CLEANUP!
#rm -Rf ${TMPDIR}


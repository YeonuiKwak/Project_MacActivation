#!/bin/bash
#Extract 3'UTR sequence inclusive in distal tandem 3'utr, not in proximal.
#Extract 3'UTR inclusive only in distal UTR.
bedtools subtract -a custom.3UTR.distal.bed12 -b custom.3UTR.proximal.bed12 -s >distal.only.UTR.region.bed
#Extract The region in fastsa.
bedtools getfasta -fi ~/Work/shared/ref/hg38/hg38.fa -bed distal.only.UTR.region.bed -split -name -s |\
      fold -w 60 >custom3UTR.distalOnly.fa

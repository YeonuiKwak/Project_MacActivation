#!/bin/bash
#Extract custom 3'UTR sequence using custom.UTR.bed12.
bedtools getfasta -fi ~/Work/shared/ref/hg38/hg38.fa -bed custom.3UTR.bed12 -split -name -s |fold -w 60 >custom3UTR.fa

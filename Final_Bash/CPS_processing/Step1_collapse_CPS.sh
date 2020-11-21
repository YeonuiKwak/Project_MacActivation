#!/bin/bash
#usage: bash step1_collapse_CPS.sh <biological rep1.bed> <biological rep2.bed> <outputName.bed>
#collapase cps positions betwee two biological replicates of CPS-seq.

cat $1 $2 \
|sort -k1,1 -k2,2n -k5,5nr |sort -k1,1 -k2,2n -u > $3

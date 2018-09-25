#!/usr/bin/env bash

fq1=sample.R1.001.fq.gz
fq2=sample.R2.001.fq.gz
OUT=sample

# takes R2 and crops last 10 bases
java -jar ${trim}/trimmomatic-0.36.jar SE -threads 20 $fq2 ${OUT}.R2.MID.fq.gz CROP:10

mid=${OUT}.R2.MID.fq.gz

# takes R2 and removes first 10 bases
java -jar $trim/trimmomatic-0.36.jar SE -threads 20 $fq2 ${OUT}.trimd.R2.fq.gz HEADCROP:10

bwa mem ${REF} ${fq1} ${OUT}.trimd.R2.fq.gz -M -t 16 -v 1 > ${OUT}.sam
#!/usr/bin/env bash

echo "file R1" > sample_R1_L001.fastq
echo "file R2" > sample_R2_L001.fastq

for f in *_R1_L001.fastq
do
fq1=$f
fq2=${f%_R1*}_R2_L001.fastq # strip shortest match of (_R1*) from f
OUT=${fq1%%_L001*} # strip longest match of (_L001*) from fq1

printf "fq1 is: %s\n" $fq1
printf "fq2 is: %s\n" $fq2
printf "out is: %s\n" $OUT

#echo $fq1
#echo $fq2
#echo $OUT
done

rm sample_R1_L001.fastq
rm sample_R2_L001.fastq

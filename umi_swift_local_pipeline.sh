#!/usr/bin/env bash

# start with a fresh session of dockers
# NOOKILLAR option for cleaning docker containers: docker rm $(docker ps -aq)

docker stop trimmomatic bwa picard snpeff gatk3.6 lofreq coveragebed
docker rm trimmomatic bwa picard snpeff gatk3.6 lofreq coveragebed

# set initial parameters
ncores=6

# TODO: vprašaj katere indel fajle se uporablja (in zakaj)

#### TOOLS ####
# TRIMMOMATIC
docker run -dt --name trimmomatic \
    -v /home/romunov/biotools/db:/db \
    -v /home/romunov/Documents/projects/swiftbio_umi_pipeline/:/reads \
    resolwebio/rnaseq:3.2.0 \
    bash --login

# BWA
docker run -dt --name bwa \
    -v /home/romunov/biotools/db/:/db \
    -v /home/romunov/Documents/projects/swiftbio_umi_pipeline/:/reads \
    resolwebio/rnaseq:3.6.1 \
    bash --login

# PICARD
docker run -dt --name picard \
    -v /home/romunov/biotools/db/:/db \
    -v /home/romunov/Documents/projects/swiftbio_umi_pipeline/:/reads \
    resolwebio/dnaseq:3.0.0 \
    bash --login

# LOFREQ
docker run -dt --name lofreq \
    -v /home/romunov/biotools/db/:/db \
    -v /home/romunov/Documents/projects/swiftbio_umi_pipeline/:/reads \
    resolwebio/dnaseq:3.1.0 \
    bash --login

# COVERAGEBED
docker run -dt --name coveragebed \
    -v /home/romunov/biotools/db/:/db \
    -v /home/romunov/Documents/projects/swiftbio_umi_pipeline/:/reads \
    resolwebio/dnaseq:3.2.1 \
    bash --login

# SNPEFF
docker run -dt --name snpeff \
    -v /home/romunov/biotools/db:/db \
    -v /home/romunov/Documents/projects/swiftbio_umi_pipeline/:/reads \
    resolwebio/legacy:1.0.0 \
    bash --login

# FGBIO (from GitHub)
fgbio=/home/romunov/biotools/fgbio/target/scala-2.12/fgbio-0.7.0-f2d5d5e-SNAPSHOT.jar

# PRIMERCLIP (latest version from GitHub)
primerclip=/home/romunov/.local/bin/primerclip

# GATK
docker run -dt --name gatk3.6 \
    -v /home/romunov/biotools/db:/db \
    -v /home/romunov/Documents/projects/swiftbio_umi_pipeline/:/reads \
    broadinstitute/genomes-in-the-cloud:2.3.1-1504795437 \
    bash --login
# usage: java -jar GATK36.jar -T RealignerTargetCreator --help

# PYTHON SCRIPT LOFREQ (from resolwe-bio package)
filterlf=/home/romunov/Documents/genialis/resolwe-bio/resolwe_bio/tools/lofreq2_indel_ovlp.py

#### TOOLS ####

#### DATABASES ####
indels_thousandg=/db/indels/1000G_phase1.indels.b37.vcf
indels_mills=/db/indels/Mills_and_1000G_gold_standard.indels.b37.vcf
# dbNSFP3.5a found here: # https://sites.google.com/site/jpopgen/dbNSFP
# database constructed using instructions here: http://snpeff.sourceforge.net/SnpSift.html#dbNSFP
dbnsfp=/db/dbNSFP/dbNSFP3.5a.txt.gz
cosmic=/db/cosmic/CosmicCodingMuts-v84.vcf
cosmicNonCod=/db/cosmic/CosmicNonCodingVariants-v84.vcf
clinvar=/db/clinvar/20180805/clinvar_20180805.vcf

REF="/db/genomes/hs_b37.fasta"
MASTER="/home/romunov/biotools/db/masterfiles/18-2132_EGFR_MID_Masterfile.txt"

# master file has been uploaded to Genialis platform for processing and bed files downloaded
BED="/db/masterfiles/18-2132_EGFR_MID_Masterfile_merged_targets_5col.bed"
BEDL="/home/romunov/biotools/db/masterfiles/18-2132_EGFR_MID_Masterfile_merged_targets_5col.bed" # bed local
NOOLAPBED="/db/masterfiles/18-2132_EGFR_MID_Masterfile_nonmerged_noolaps_targets_5col.bed"

known_dbsnp="/db/known_sites/dbsnp_138.b37.vcf"
known_1000g="/db/known_sites/1000G_phase1.snps.high_confidence.b37.vcf"
known_hapmap="/db/known_sites/hapmap_3.3.b37.vcf"
#### DATABASES ####

set -e # exit if simple command exits with non-zero status unless in a loop
set -x # print trace of sample commands and their arguments

for f in *_R1_001.fastq.gz
do

fq1=${f}
fq2=${f%_R1*}_R2_001.fastq.gz # strips string of _R1*
OUT=${fq1%%_L001*} # strips string of _L001*

echo $(printf "##### PROCESSING SAMPLE %s" ${OUT})

########## PRE-MID analysis & POST-MID with fgbio ##############
# CROP: Cut the read to a specified length by removing bases from the end - produce a MID file
run_trimmomatic1=$(printf "TrimmomaticSE -threads %s %s %s.R2.MID.fq.gz CROP:10" ${ncores} ${fq2} ${OUT})
time docker exec -w /reads -t trimmomatic bash --login  $run_trimmomatic1

# MID fastq file
mid=${OUT}.R2.MID.fq.gz

# HEADCROP: Cut the specified number of bases from the start of the read
run_trimmomatic2=$(printf "TrimmomaticSE -threads %s %s %s.trimd.R2.fq.gz HEADCROP:10" ${ncores} ${fq2} ${OUT})
time docker exec -w /reads -t trimmomatic bash --login $run_trimmomatic2

# align using bwa with verbosity 1 (to print errors only), mark sec alignments and add RG
#bwa mem ${REF} ${fq1} ${OUT}.trimd.R2.fq.gz -M -t 16 -v 1 > ${OUT}.sam # TODO: why aligning raw R1 and trimmed R2 reads?
run_bwa1=$(printf "bwa mem %s %s %s.trimd.R2.fq.gz -M -t %s -v 1 > %s.sam" ${REF} ${fq1} ${OUT} ${ncores} ${OUT})
time docker exec -w /reads -t bwa bash --login -c "$run_bwa1"

# Sort, Add RG for nopclip bam
run_picard_addreplace_readgrp=$(printf \
"picard-tools AddOrReplaceReadGroups \
    I=%s.sam \
    O=%s.preMID.bam \
    SO=coordinate \
    VALIDATION_STRINGENCY=STRICT \
    RGID=Accel \
    RGLB=MID \
    RGSM=%s \
    RGPL=Illumina \
    RGPU=Miseq \
    CREATE_INDEX=TRUE" \
${OUT} ${OUT} ${OUT})

time docker exec -w /reads -t picard bash --login $run_picard_addreplace_readgrp

# query sort for primerclip
#java -jar $picard SortSam I=$OUT.sam O=$OUT.qsort.sam SO=queryname
run_sortsam1=$(printf \
    "picard-tools SortSam \
    I=%s.sam \
    O=%s.qsort.sam \
    SO=queryname" ${OUT} ${OUT})

time docker exec -w /reads -t picard bash --login $run_sortsam1

#PrimerCLIP before MID
#/home/sandhu/.local/bin/primerclip $MASTER $OUT.qsort.sam preM${OUT}.preMID.pclip.sam
#/home/smrtanalysis/.local/bin/primerclip $MASTER $OUT.qsort.sam preM${OUT}.preMID.pclip.sam

# Sort
#java -jar $picard AddOrReplaceReadGroups I=preM${OUT}.preMID.pclip.sam O=$OUT.preMID.pclip.bam \
# SO=coordinate VALIDATION_STRINGENCY=STRICT RGID=Accel RGLB=MID RGSM=$OUT RGPL=Illumina \
# RGPU=Miseq CREATE_INDEX=TRUE

## variant call preMID primer trimmed bam
#lofreq call --call-indels -f $REF $OUT.preMID.bam -l $BED -o $OUT.preMID.noPclip.lf.vcf
run_lf1=$(printf \
"lofreq call-parallel \
    --pp-threads %s \
    --call-indels \
    --ref %s \
    %s.preMID.bam \
    --bed %s \
    --out %s.preMID.noPclip.lf.vcf" \
    ${ncores} \
    ${REF} \
    ${OUT} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t lofreq bash --login -c "$run_lf1"

### MID processing ########
# annotate bam with MIDs from the I2 file
time java -Xmx8g -jar $fgbio AnnotateBamWithUmis -i $OUT.preMID.bam -f $mid -o $OUT.bamAnnotdWumi.bam

#bedtools to get per amplicon coverage
# NOTE: use of parameters in coverageBed has changed with v2.24. Use -abam for the bed file and -b for .bam file.
#time coverageBed -abam $OUT.bamAnnotdWumi.bam -b ${BEDL} > $OUT.bamAnnotdWumi.cov
#time coverageBed -abam $OUT.preMID.bam -b ${BEDL} -d > $OUT.preMID.covd
#time coverageBed -abam $OUT.preMID.bam -b ${BEDL} > $OUT.preMID.cov

run_cb1=$(printf "coverageBed -b %s.bamAnnotdWumi.bam -abam %s > %s.bamAnnotdWumi.cov" ${OUT} ${BED} ${OUT})
run_cb2=$(printf "coverageBed -b %s.preMID.bam -abam %s -d > %s.preMID.covd" ${OUT} ${BED} ${OUT})
run_cb3=$(printf "coverageBed -b %s.preMID.bam -abam %s > %s.preMID.cov" ${OUT} ${BED} ${OUT})

time docker exec -w /reads -t coveragebed bash --login -c "$run_cb1"
time docker exec -w /reads -t coveragebed bash --login -c "$run_cb2"
time docker exec -w /reads -t coveragebed bash --login -c "$run_cb3"

#RevertSam to sanitise
#java -jar $picard RevertSam I=$OUT.bamAnnotdWumi.bam O=$OUT.sanitised.bam \
#    SANITIZE=true REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false
run_revertsam=$(printf \
"picard-tools RevertSam \
    I=%s.bamAnnotdWumi.bam \
    O=%s.sanitised.bam \
    SANITIZE=true \
    REMOVE_DUPLICATE_INFORMATION=false \
    REMOVE_ALIGNMENT_INFORMATION=false" \
    ${OUT} ${OUT})
time docker exec -w /reads -t coveragebed bash --login -c "$run_revertsam"

#SetMateInfo to bring MQ tags
time java -jar $fgbio SetMateInformation -i $OUT.sanitised.bam -o $OUT.sanitised.setmateinfo.bam

#sort by queryname
#java -jar $picard SortSam I=$OUT.sanitised.setmateinfo.bam O=$OUT.sanitised.setmateinfo.qsort.bam SO=queryname
run_sortsam2=$(printf \
"picard-tools SortSam \
    I=%s.sanitised.setmateinfo.bam \
    O=%s.sanitised.setmateinfo.qsort.bam \
    SO=queryname"\
     ${OUT} ${OUT})

time docker exec -w /reads -t picard bash --login -c "$run_sortsam2"

# Group Reads by UMI
time java -jar $fgbio GroupReadsByUmi -s adjacency --edits 1 -i $OUT.sanitised.setmateinfo.qsort.bam -o $OUT.groupdbyMID.bam

# Make consensus N molecules to make consensus Out is coordinate sorted
time java -jar $fgbio CallMolecularConsensusReads -M 3 -i $OUT.groupdbyMID.bam -o $OUT.M3.consensus.bam -r $OUT.M3.notused4consensus.bam

# argument min reads required, so I use the "default" 1
time java -jar $fgbio CallMolecularConsensusReads -M 1 -i $OUT.groupdbyMID.bam -o $OUT.consensus.bam -r $OUT.notused4consensus.bam

# Extract paired end reads from bam to fastq files
#java -Xmx128g -jar $picard SamToFastq I=$OUT.M3.consensus.bam F=$OUT.consensus.R1.fq F2=$OUT.consensus.R2.fq FU=$OUT.unpaired.fq
run_samtofq=$(printf \
"picard-tools SamToFastq \
    I=%s.M3.consensus.bam \
    F=%s.consensus.R1.fq \
    F2=%s.consensus.R2.fq \
    FU=%s.unpaired.fq" \
    ${OUT} \
    ${OUT} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads picard bash --login -c "$run_samtofq"

# Realign
#bwa mem $REF $OUT.consensus.R1.fq $OUT.consensus.R2.fq -t 24 -M > ${OUT}_consensus.sam
##bwa mem $REF $OUT.consensus.R1.fq $OUT.consensus.R2.fq -t 24 -M -I 100,50,50,500 > ${OUT}_consensus.sam
##bwa mem $REF $OUT.pclip.consensus.R1.fq $OUT.pclip.consensus.R2.fq -t 24 -M > ${OUT}_pclip.consensus.sam
run_realign=$(printf \
"bwa mem \
    %s \
    %s.consensus.R1.fq \
    %s.consensus.R2.fq \
    -t %s \
    -M > %s_consensus.sam" \
    ${REF} \
    ${OUT} \
    ${OUT} \
    ${ncores} \
    ${OUT})

time docker exec -w /reads -t bwa bash --login -c "$run_realign"

## qsort for primer-clip
#java -jar $picard SortSam I=${OUT}_consensus.sam O=$OUT.consensus.qsort.sam SO=queryname
run_sortsam3=$(printf \
"picard-tools SortSam \
    I=%s_consensus.sam \
    O=%s.consensus.qsort.sam \
    SO=queryname" \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t picard bash --login -c "$run_sortsam3"

#PrimerCLIP
time primerclip $MASTER $OUT.consensus.qsort.sam ${OUT}.fgbio.pclip.sam
#primerclip $MASTER $OUT.pclip.consensus.qsort.sam pclip${OUT}.pclip.fgbio.pclip.sam

#sort and Add read groups
#java -Xmx64g -jar $picard AddOrReplaceReadGroups I=${OUT}.fgbio.pclip.sam \
#    O=$OUT.fgbio.bam RGID=MID-Amp RGLB=$OUT RGSM=NA12878 RGPL=Illumina \
#    RGPU=MiSeq SO=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=STRICT
# TODO: how to define RGSM value? based on description from the manual, it feels like this is something arbitrary
run_picard_oddreplace_readgrp2=$(printf \
"picard-tools AddOrReplaceReadGroups \
    I=%s.fgbio.pclip.sam \
    O=%s.fgbio.bam \
    RGID=MID-Amp \
    RGLB=%s \
    RGSM=NA12878 \
    RGPL=Illumina \
    RGPU=MiSeq \
    SO=coordinate \
    CREATE_INDEX=TRUE \
    VALIDATION_STRINGENCY=STRICT" \
    ${OUT} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t picard bash --login -c "$run_picard_oddreplace_readgrp2"

# make BED file for picard
samtools view -H $OUT.fgbio.bam > header.txt
cat header.txt /home/romunov/biotools$BED > $PWD/BED.picard.bed
pBED=BED.picard.bed

# Collect TargetedPCRMetrics from the TWO BAMs #######
#java -Xmx65g -jar $picard CollectTargetedPcrMetrics I=$OUT.preMID.bam \
#    O=$OUT.preMID.targetPCRmetrics.txt TI=$pBED AI=$pBED R=$REF \
#    PER_TARGET_COVERAGE=$OUT.preMID.perTargetCov.txt
run_coltarmet1=$(printf \
"picard-tools CollectTargetedPcrMetrics \
    I=%s.preMID.bam \
    O=%s.preMID.targetPCRmetrics.txt \
    TI=%s \
    AI=%s \
    R=%s \
    PER_TARGET_COVERAGE=%s.preMID.perTargetCov.txt" \
    ${OUT} \
    ${OUT} \
    ${pBED} \
    ${pBED} \
    ${REF} \
    ${OUT})

time docker exec -w /reads -t picard bash --login -c "$run_coltarmet1"

#java -Xmx64g -jar $picard CollectTargetedPcrMetrics I=$OUT.fgbio.bam \
#    O=$OUT.fgbio.targetPCRmetrics.txt TI=$pBED AI=$pBED R=$REF \
#    PER_TARGET_COVERAGE=$OUT.fgbio.perTargetCov.txt
run_coltarmet2=$(printf \
"picard-tools CollectTargetedPcrMetrics \
    I=%s.fgbio.bam \
    O=%s.fgbio.targetPCRmetrics.txt \
    TI=%s \
    AI=%s \
    R=%s \
    PER_TARGET_COVERAGE=%s.fgbio.perTargetCov.txt" \
    ${OUT} \
    ${OUT} \
    ${pBED} \
    ${pBED} \
    ${REF} \
    ${OUT})

time docker exec -w /reads -t picard bash --login -c "$run_coltarmet2"

#BEDtools
#coverageBed -abam $OUT.fgbio.bam -b $BED > $OUT.fgbio.cov
#coverageBed -abam $OUT.fgbio.bam -b $BED -d > $OUT.fgbio.covd
run_cb4=$(printf "coverageBed -b %s.fgbio.bam -abam %s > %s.fgbio.cov" ${OUT} ${BED} ${OUT})
run_cb5=$(printf "coverageBed -b %s.fgbio.bam -abam %s -d > %s.fgbio.covd" ${OUT} ${BED} ${OUT})

#coverageBed -abam $OUT.fgbio.bam -b $NOOLAPBED > $OUT.noolap.fgbio.cov
#coverageBed -abam $OUT.fgbio.bam -b $NOOLAPBED -d > $OUT.noolap.fgbio.covd
run_cb6=$(printf "coverageBed -b %s.fgbio.bam -abam %s > %s.noolap.fgbio.cov" ${OUT} ${NOOLAPBED} ${OUT})
run_cb7=$(printf "coverageBed -b %s.fgbio.bam -abam %s -d > %s.noolap.fgbio.covd" ${OUT} ${NOOLAPBED} ${OUT})

# pre MID cov covd files with noolap BED
#coverageBed -abam $OUT.preMID.bam -b $NOOLAPBED > $OUT.noolap.preMID.cov
#coverageBed -abam $OUT.preMID.bam -b $NOOLAPBED -d > $OUT.noolap.preMID.covd
run_cb8=$(printf "coverageBed -b %s.preMID.bam -abam %s > %s.noolap.preMID.cov" ${OUT} ${NOOLAPBED} ${OUT})
run_cb9=$(printf "coverageBed -b %s.preMID.bam -abam %s -d > %s.noolap.preMID.covd" ${OUT} ${NOOLAPBED} ${OUT})

time docker exec -w /reads -t coveragebed bash --login -c "$run_cb4"
time docker exec -w /reads -t coveragebed bash --login -c "$run_cb5"
time docker exec -w /reads -t coveragebed bash --login -c "$run_cb6"
time docker exec -w /reads -t coveragebed bash --login -c "$run_cb7"
time docker exec -w /reads -t coveragebed bash --login -c "$run_cb8"
time docker exec -w /reads -t coveragebed bash --login -c "$run_cb9"


############# BQSR fgbio bam      ##########
#Create target interval for Indelrealigner
#$GATK  -T RealignerTargetCreator -R $REF \
#    -known $indel_thousandg \
#    -known $indel_mills \
#    -I $OUT.fgbio.bam -L $BED -o ${OUT}.forIndelRealigner.intervals

# in case index is missing, create one
#docker exec -w /reads -t picard bash --login picard-tools CreateSequenceDictionary R=/db/genomes/hs_b37.fasta O=/db/genomes/hs_b37.dict
run_gatk1=$(printf \
"java -jar /usr/gitc/GATK36.jar -T RealignerTargetCreator \
    -R %s \
    -known %s \
    -known %s \
    -I %s.fgbio.bam \
    -L %s \
    -o %s.forIndelRealigner.intervals" \
    ${REF} \
    ${indels_thousandg} \
    ${indels_mills} \
    ${OUT} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk1"

#GATK Indel Realignment
#$GATK  -T IndelRealigner -R $REF -I $OUT.fgbio.bam \
#    -known $indels_thousandg \
#    -known $indels_mills \
#    --targetIntervals ${OUT}.forIndelRealigner.intervals \
#    -L $BED -o ${OUT}.realigned.bam
run_gatk2=$(printf \
"java -jar /usr/gitc/GATK36.jar -T IndelRealigner \
    -R %s\
    -I %s.fgbio.bam \
    -known %s \
    -known %s \
    --targetIntervals %s.forIndelRealigner.intervals \
    -L %s \
    -o %s.realigned.bam" \
    ${REF} \
    ${OUT} \
    ${indels_thousandg} \
    ${indels_mills} \
    ${OUT} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk2"

#$GATK  -T BaseRecalibrator -R $REF -I $OUT.fgbio.bam \
#    -nct 12 \
#    --knownSites $anno/data/GRCh37.75/dbsnp141.vcf \
#    --knownSites $anno/data/GRCh37.75/1000Genome_extra.vcf \
#    --knownSites $anno/data/GRCh37.75/1000G_phase1.snps.high_confidence.b37.vcf \
#    --knownSites $anno/data/GRCh37.75/hapmap_3.3.b37.vcf \
#    -L $BED -o ${OUT}.recal_data.table
# Recalibrate base quality scores.
run_gatk3=$(printf \
"java -jar /usr/gitc/GATK36.jar -T BaseRecalibrator \
    -R %s\
    -I %s.fgbio.bam \
    -nct %s \
    --knownSites %s \
    --knownSites %s \
    --knownSites %s \
    -L %s \
    -o %s.recal_data.table" \
    ${REF} \
    ${OUT} \
    ${ncores} \
    ${known_dbsnp} \
    ${known_1000g} \
    ${known_hapmap} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk3"

#generate Recalibrated bam
#$GATK  -T PrintReads -R $REF -I ${OUT}.realigned.bam \
#    -BQSR ${OUT}.recal_data.table -o ${OUT}.fgbio.bqsrCal.bam
# Create recalibrated BAM file.
run_gatk4=$(printf \
"java -jar /usr/gitc/GATK36.jar  -T PrintReads \
    -R %s \
    -I %s.realigned.bam \
    --BQSR %s.recal_data.table \
    -o %s.fgbio.bqsrCal.bam" \
    ${REF} \
    ${OUT} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk4"

#########################################
#### BQSR pre_MID bam####################
#Create target interval for Indelrealigner
#$GATK  -T RealignerTargetCreator -R $REF \
#    -known $indels_thousandg \
#    -known $indels_mills \
#    -I $OUT.preMID.bam -L $BED -o ${OUT}.forIndelRealigner.intervals
run_gatk5=$(printf \
"java -jar /usr/gitc/GATK36.jar \
    -T RealignerTargetCreator \
    -R %s \
    -known %s \
    -known %s\
    -I %s.preMID.bam \
    -L %s \
    -o %s.forIndelRealigner.intervals" \
    ${REF} \
    ${indels_thousandg} \
    ${indels_mills} \
    ${OUT} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk5"

#GATK Indel Realignment
#$GATK  -T IndelRealigner -R $REF -I $OUT.preMID.bam \
#    -known $indels_thousandg \
#    -known $indels_mills \
#    --targetIntervals ${OUT}.forIndelRealigner.intervals \
#    -L $BED -o ${OUT}.preMID.realigned.bam
# Realign reads such that the number of mismatching bases is minimized across all reads.
run_gatk6=$(printf \
"java -jar /usr/gitc/GATK36.jar \
    -T IndelRealigner \
    -R %s \
    -I %s.preMID.bam \
    -known %s\
    -known %s \
    --targetIntervals %s.forIndelRealigner.intervals \
    -L %s \
    -o %s.preMID.realigned.bam" \
    ${REF} \
    ${OUT} \
    ${indels_thousandg} \
    ${indels_mills} \
    ${OUT} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk6"

#$GATK  -T BaseRecalibrator -R $REF -I $OUT.preMID.bam \
#    --knownSites $anno/data/GRCh37.75/dbsnp141.vcf -nct 12 \
#    --knownSites $anno/data/GRCh37.75/1000Genome_extra.vcf \
#    --knownSites $anno/data/GRCh37.75/1000G_phase1.snps.high_confidence.b37.vcf \
#    --knownSites $anno/data/GRCh37.75/hapmap_3.3.b37.vcf \
#    -L $BED -o ${OUT}.recal_data.table

time run_gatk7=$(printf \
"java -jar /usr/gitc/GATK36.jar \
    -T BaseRecalibrator \
    -R %s \
    -I %s.preMID.bam \
    --knownSites %s \
    --knownSites %s \
    --knownSites %s \
    -L %s \
    -o %s.recal_data.table" \
    ${REF} \
    ${OUT} \
    ${known_dbsnp} \
    ${known_1000g} \
    ${known_hapmap} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk7"

#generate Recalibrated bam
#$GATK  -T PrintReads -R $REF -I ${OUT}.preMID.realigned.bam -BQSR ${OUT}.recal_data.table -o ${OUT}.preMID.bqsrCal.bam
run_gatk8=$(printf \
"java -jar /usr/gitc/GATK36.jar \
    -T PrintReads \
    -R %s \
    -I %s.preMID.realigned.bam \
    --BQSR %s.recal_data.table \
    -o %s.preMID.bqsrCal.bam" \
    ${REF} \
    ${OUT} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk8"

################################################################
###############################################################

#low freq variant calls
#lofreq call --call-indels -f $REF ${OUT}.fgbio.bqsrCal.bam -l $BED -o $OUT.fgbio.bqsr.lf.vcf
run_lf2=$(printf \
"lofreq call-parallel \
    --pp-threads %s \
    --call-indels \
    --ref %s \
    %s.fgbio.bqsrCal.bam \
    --bed %s \
    --out %s.fgbio.bqsr.lf.vcf"\
    ${ncores} \
    ${REF} \
    ${OUT} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t lofreq bash --login -c "$run_lf2"

#lofreq call --call-indels -f $REF ${OUT}.preMID.bqsrCal.bam -l $BED -o $OUT.preMID.bqsr.lf.vcf
run_lf3=$(printf \
"lofreq call-parallel \
    --pp-threads %s \
    --call-indels \
    --ref %s \
    %s.preMID.bqsrCal.bam \
    --bed %s \
    --out %s.preMID.bqsr.lf.vcf" \
    ${ncores} \
    ${REF} \
    ${OUT} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t lofreq bash --login -c "$run_lf3"

#lofreq call --call-indels -f $REF ${OUT}.fgbio.bam -l $BED -o $OUT.fgbio.NObqsr.lf.vcf
run_lf4=$(printf \
"lofreq call-parallel \
    --pp-threads %s \
    --call-indels \
    --ref %s \
    %s.fgbio.bam \
    --bed %s \
    --out %s.fgbio.NObqsr.lf.vcf" \
    ${ncores} \
    ${REF} \
    ${OUT} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t lofreq bash --login -c "$run_lf4"

#lofreq vcf merge overlapping indels
time python ${filterlf} ${OUT}.fgbio.bqsr.lf.vcf ${OUT}.fgbio.bqsr.lf.mergeIndels.vcf
time python ${filterlf} ${OUT}.preMID.bqsr.lf.vcf ${OUT}.preMID.bqsr.lf.mergeIndels.vcf

#germline var calling GATK HC
#$GATK -T HaplotypeCaller -I $OUT.fgbio.bqsrCal.bam -R $REF -L $BED -o $OUT.fgbio.bqsr.gatkHC.vcf
run_gatk9=$(printf \
"java -jar /usr/gitc/GATK36.jar \
    -T HaplotypeCaller \
    -I %s.fgbio.bqsrCal.bam \
    -R %s \
    -L %s \
    -o %s.fgbio.bqsr.gatkHC.vcf" \
    ${OUT} \
    ${REF} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk9"

#$GATK -T HaplotypeCaller -I $OUT.preMID.bqsrCal.bam -R $REF -L $BED -o $OUT.preMID.bqsr.gatkHC.vcf
run_gatk10=$(printf \
"java -jar /usr/gitc/GATK36.jar \
    -T HaplotypeCaller \
    -I %s.preMID.bqsrCal.bam \
    -R %s \
    -L %s \
    -o %s.preMID.bqsr.gatkHC.vcf" \
    ${OUT} \
    ${REF} \
    ${BED} \
    ${OUT})

time docker exec -w /reads -t gatk3.6 bash --login -c "$run_gatk10"

#java -Xmx64g -jar $anno/SnpSift.jar annotate $cosmic $OUT.fgbio.bqsr.lf.mergeIndels.vcf > $OUT.fgbio.lf.cosmic.vcf
run_snpeff1=$(printf \
"SnpSift \
    annotate %s \
    %s.fgbio.bqsr.lf.mergeIndels.vcf > %s.fgbio.lf.cosmic.vcf" \
    ${cosmic} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t snpeff bash --login -c "$run_snpeff1"

#java -Xmx64g -jar $anno/SnpSift.jar annotate $cosmic $OUT.fgbio.bqsr.gatkHC.vcf > $OUT.fgbio.gatkHC.cosmic.vcf
run_snpeff2=$(printf \
"SnpSift \
    annotate %s \
    %s.fgbio.bqsr.gatkHC.vcf > %s.fgbio.gatkHC.cosmic.vcf" \
    ${cosmic} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t snpeff bash --login -c "$run_snpeff2"

##java -Xmx64g -jar $anno/SnpSift.jar annotate $clinvar $OUT.fgbio.lf.cosmic.vcf > $OUT.fgbio.lf.cosmic.clinvar.vcf
run_snpeff3=$(printf \
"SnpSift \
    annotate %s \
    %s.fgbio.lf.cosmic.vcf > %s.fgbio.lf.cosmic.clinvar.vcf" \
    ${clinvar} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t snpeff bash --login -c "$run_snpeff3"

##java -Xmx64g -jar $anno/SnpSift.jar annotate $clinvar $OUT.fgbio.gatkHC.cosmic.vcf > $OUT.fgbio.gatkHC.cosmic.clinvar.vcf
run_snpeff4=$(printf \
"SnpSift \
    annotate %s \
    %s.fgbio.gatkHC.cosmic.vcf > %s.fgbio.gatkHC.cosmic.clinvar.vcf" \
    ${clinvar} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t snpeff bash --login -c "$run_snpeff4"

# dbNSFP annotation for polyphen and pathogenic score
#java -Xmx64g -jar $anno/SnpSift.jar dbNSFP -db $dbnsfp $OUT.fgbio.lf.cosmic.clinvar.vcf > $OUT.fgbio.lf.cosmic.clinvar.dbnsfp.vcf
run_snpeff5=$(printf \
"SnpSift \
    dbnsfp \
    -db %s \
    %s.fgbio.lf.cosmic.clinvar.vcf > %s.fgbio.lf.cosmic.clinvar.dbnsfp.vcf" \
    ${dbnsfp} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t snpeff bash --login -c "$run_snpeff5"

#java -Xmx64g -jar $anno/SnpSift.jar dbNSFP -db $dbnsfp $OUT.fgbio.gatkHC.cosmic.clinvar.vcf > $OUT.fgbio.gatkHC.cosmic.clinvar.dbnsfp.vcf
run_snpeff6=$(printf \
"SnpSift \
    dbnsfp \
    -db %s \
    %s.fgbio.gatkHC.cosmic.clinvar.vcf > %s.fgbio.gatkHC.cosmic.clinvar.dbnsfp.vcf" \
    ${dbnsfp} \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t snpeff bash --login -c "$run_snpeff6"

# Final variants from clinvar, cosmic &dbnsfp
#java -Xmx64g -jar $anno/SnpSift.jar extractFields $OUT.fgbio.lf.cosmic.clinvar.dbnsfp.vcf \
# CHROM POS ID REF ALT QUAL DP AF SB DP4 "CLNDN" "CLNSIG" "CLNVI" "CLNVC" "GENEINFO" "ORIGIN" "RS" "MC" "dbNSFP_LRT_pred" \
#    "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_MutationTaster_pred" "AF_ESP" "AF_TGP" > $OUT.fgbio.lf.finalvars.txt
run_final1=$(printf \
"SnpSift \
    extractFields \
    %s.fgbio.lf.cosmic.clinvar.dbnsfp.vcf \
    CHROM POS ID REF ALT QUAL DP AF SB DP4 \"CLNDN\" \"CLNSIG\" \"CLNVI\" \"CLNVC\" \"GENEINFO\" \"ORIGIN\" \"RS\" \"MC\" \
    \"dbNSFP_LRT_pred\" \"dbNSFP_Polyphen2_HDIV_pred\" \"dbNSFP_MutationTaster_pred\" \"AF_ESP\" \"AF_TGP\" \
    > %s.fgbio.lf.finalvars.txt" \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t snpeff bash --login -c "$run_final1"

#time java -Xmx64g -jar $anno/SnpSift.jar extractFields $OUT.fgbio.gatkHC.cosmic.clinvar.dbnsfp.vcf \
# CHROM POS ID REF ALT QUAL DP AF SB DP4 "CLNDN" "CLNSIG" "CLNVI" "CLNVC" "GENEINFO" "ORIGIN" "RS" "MC" "dbNSFP_LRT_pred" \
#    "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_MutationTaster_pred" "AF_ESP" "AF_TGP" > $OUT.fgbio.gatkHC.finalvars.txt
run_final2=$(printf \
"SnpSift \
    extractFields \
    %s.fgbio.gatkHC.cosmic.clinvar.dbnsfp.vcf \
    CHROM POS ID REF ALT QUAL DP AF SB \"CLNDN\" \"CLNSIG\" \"CLNVI\" \"CLNVC\" \"GENEINFO\" \"ORIGIN\" \"RS\" \"MC\" \
    \"dbNSFP_LRT_pred\" \"dbNSFP_Polyphen2_HDIV_pred\" \"dbNSFP_MutationTaster_pred\" \"AF_ESP\" \"AF_TGP\" \
    > %s.fgbio.gatkHC.finalvars.txt" \
    ${OUT} \
    ${OUT})

time docker exec -w /reads -t snpeff bash --login -c "$run_final2"

done

######  clean up ######
rm *.MID.fq.gz
#rm unpaired*gz
rm *var.vcf
rm *cosmic.vcf
#rm *sus.bam
#rm *MID.bam
rm *.sanitise*.bam
#rm *notused4consensus.bam
#rm *.realigned.ba[mi]
rm *table

#############################################################
##########   summarize coverage results and reporting ########
for f in *covd
do
    # TODO: m*0.2 ?????
    echo $(printf "Processing %s" ${f})
    awk '{sum+=$7}END{m=(sum/NR); b=m*0.2; print m, b}' ${f} > ${f}_covd.tmp
    awk 'BEGIN{n=0}NR==FNR{m=$1;b=$2;next}{if($7>=b)n++}END{print m,(n/FNR*100.0)}' \
     OFS="\t" ${f}_covd.tmp ${f} > ${f%.covd}_covMetrics.txt
done

rm *covd.tmp

#awk 'BEGIN{n=0}NR==FNR{m=$1;b=$2;next}{if($6>=b)n++}END{print m,b,(n/FNR*100.0)}'
# Summarize PCR metrics
# new PICARD rs2

for f in *targetPCRmetrics.txt
do
    echo $(printf "Processing %s" ${f})
    awk -v n=${f%%.target*} 'NR==1{print n,$5,$11,$14*100,$22*100.0,($17+$18)/$7*100.0}' \
        OFS="\t" <(head -n8 $f | tail -n1) > ${f%%.txt}_summary.txt
    f2=${f%%.target*}_covMetrics.txt
    paste ${f%%.txt}_summary.txt $f2 > ${f2%%_cov*}_combined_cov_metrics.txt
done

echo "SampleID    Total_Reads    #UQ_Reads_Aligned    %UQ_Reads_Aligned    %Bases_OnTarget_Aligned     %Bases_OnTarget_Total    Mean_Coverage    %Coverage_Uniformity" \
    > final_metrics_report.txt
cat *_combined_cov_metrics.txt >> final_metrics_report.txt

rm *_summary.txt
rm *_covMetrics.txt
rm *cov_metrics.txt
rm masterparsefails.log
rm BED.picard.bed

# run quality check statistics on outputs of all samples
./mid_QC_module.sh

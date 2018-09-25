#! /bin/bash
# S Sandhu 180725
# OUTPUTs for MID amplicon QC: tables and tables to be plotted

###########################################################################################################
######### pre and post MID read coverage TABLE to ensure sufficient seq depth was obtained ################
## Table 1: Per target Read coverage with and without MIDs ##########
## help text: To ensure sufficient seq depth was obtained.?!###

printf "Chr\tStart\tEND\tAmpliconID\tLength\n" > rowheader

#make the header column file
for f in *olap.preMID.cov
do
awk '{print $1,$2,$3,$5,$7}' OFS="\t" $f > colheader
done 

cat rowheader colheader > cov_table_header

## parse read coverage pre and post MID from bedtools output 
for f in *olap.preMID.cov
do
#g=${f%.no*}.noolap.fgbio.cov
echo "${f%.no*}_preMID" > ${f}.preCov.tmp
awk '{print $6}' $f >> ${f}.preCov.tmp

echo "${f%.no*}_postMID" > ${f}.postCov.tmp
awk '{print $6}' ${f%.no*}.noolap.fgbio.cov >> ${f}.postCov.tmp

paste ${f}.preCov.tmp ${f}.postCov.tmp > ${f}.pre-postCov

done

paste cov_table_header *pre-postCov > pre_post-MID_coverage.txt

############################################################################
######### MID family size distribution i.e. #duplicates per MID #############
## Table 2: MID family size distribution (Overall Average fam size and Real family size (at least 3 molecules per MID)

for f in *groupdbyMID.bam
do
samtools view $f |grep -o "RX:Z.*." | sed 's/RX:Z://g' | sort | uniq -c | sort -k1,1n | sed -e 's/^\s*//g' -e 's/\s/\t/g' > ${f%.group*}.mid_fam_counts.txt

## *2hist file needs plotted into a histogram (may be some Rscript or python)
cut -f1 ${f%.group*}.mid_fam_counts.txt > ${f%.group*}.counts2hist
cat ${f%.group*}.counts2hist | sort | uniq -c |sort -k2,2n | sed -e 's/^\s*//g' -e 's/\s/\t/g' > ${f%.group*}.freq.hist.txt

#printf "Average (overall) Family Size\n"
echo ${f%.group*} > ${f%.group*}.avgFamSize
awk '{sum+=$1}END{print sum/NR}' ${f%.group*}.mid_fam_counts.txt >> ${f%.group*}.avgFamSize

#printf "Average (>2 members) Family Size\n"
awk '{if($1>3) print}' ${f%.group*}.mid_fam_counts.txt| awk '{sum+=$1}END{print sum/NR}' > ${f%.group*}.realFamSize
cat ${f%.group*}.avgFamSize ${f%.group*}.realFamSize > ${f%.group*}.avg_real_famSize.txt

done

printf "Sample\nAverage (overall) Family Size\nAverage (>2 members) Family Size\n" > fam_size_header
paste fam_size_header *famSize.txt > MID_FamilySize_Stats_report.txt

## Figure 1: Histogram of MID Family size distribution
Rscript ./R/plot_family_size_histogram.R

rm *tmp
rm *counts2hist

###############################################################################################
###########   Variant binning in various AF buckets pre-postMID   #############################
### Classify variants into 6 Allele frequency bins from pre & post MID vcf #############
###  to show noise clean up due to false positive reduction when using MIDs #################
# Table 3: False positive variant clean up with MIDs

## use bqsr lofreq mergeindels vcf
for f in *.bqsr.lf.mergeIndels.vcf
do

## total variants preMID
echo ${f%.bqsr*} > ${f%.bqsr}.totalvars.tmp
grep -vc ^"#" $f >> ${f%.bqsr}.totalvars.tmp

## total variants postMID/fgbio
#echo ${f%.bqsr*} > ${f%.bqsr}.totalvars.tmp
#grep -vc ^"#" $f >> ${f%.bqsr}.totalvars.tmp

## bin variants in AF <0.1 
grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if($9<0.00099) n++}END {print n}' > ${f%.bqsr}.af0_lt0.1p.tmp

## bin variants in AF < 0.5%
grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.00099) && ($9<0.0049)) n++}END {print n}' > ${f%.bqsr}.af1_0.1-.5p.tmp

### bin variants with AF between 0.5 and 1%
grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.0049)&&($9<0.0099)) n++}END {print n}' > ${f%.bqsr}.af2_0.5-1p.tmp

### bin variants with AF between 1% and 5%
grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.0099)&&($9<0.0499)) n++}END {print n}' > ${f%.bqsr}.af3_1-5p.tmp

### bin variants with AF between 5% and 10%
grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.0499)&&($9<0.0999)) n++}END {print n}' > ${f%.bqsr}.af4_5-10p.tmp

### bin variants with AF between 10% and 50%
grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if(($9>0.0999)&&($9<0.4999)) n++}END {print n}' > ${f%.bqsr}.af5_10-50p.tmp

### bin variants with AF >50%
grep -v ^"#" $f | sed -e 's/;/\t/g' -e 's/AF=//g'| awk 'BEGIN{n=0} {if($9>=0.4999) n++}END {print n}' > ${f%.bqsr}.af6_gt50p.tmp

cat ${f%.bqsr}.totalvars.tmp ${f%.bqsr}.af*tmp > ${f%.bqsr}_varbins
done

#### merge from all samples for final report/Table
printf "Sample\nTotal_Vars\n<0.1%%\n0.1-0.5%%\n0.5-1%%\n1-5%%\n5-10%%\n10-50%%\n>= 50%%\n" > bin_header

paste bin_header *_varbins > variant_AFbins_prePOST-MID.txt

## Figure 2: Variants in allele frequencies
Rscript ./R/plot_variant_af_bins.R --input variant_AFbins_prePOST-MID.txt

#clean up
rm *tmp
rm *bins
rm *header
rm *Size
rm *Size.txt
rm *counts.txt
rm *pre-postCov

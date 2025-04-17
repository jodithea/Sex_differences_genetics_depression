#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=20GB
#PBS -l ncpus=1
#PBS -J 0-1

# This script formats input GWAS sumstats
# need to be in the GCTA-COJO .ma format = header row of SNP A1 A2 freq b se p N (SNP identifier (rsID), the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size)


### Environment ###




### Preamble ###

directory=/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/

input=(/path/UKBiobank_combined2/sex_stratified_downsampled/GWAS_UKB_depression_females_autosomes_and_X_QCed_INFO0.6_MAF0.01.txt /path/UKBiobank_combined2/sex_stratified_downsampled/GWAS_UKB_depression_males_autosomes_and_X_QCed_INFO0.6_MAF0.01.txt)

name=$(echo ${input[${PBS_ARRAY_INDEX}]} | sed 's|/path/UKBiobank_combined2/sex_stratified_downsampled/||' | sed 's|_QCed_INFO0.6_MAF0.01.txt||')



### Run script ###

cd ${directory}


awk 'NR==1 {print "SNP A1 A2 freq b se p N"} \
NR>1 {print $1, $4, $5, $6, $7, $8, $9, $10}' \
${input[${PBS_ARRAY_INDEX}]} \
> ${name}_formatted_SBayes.ma



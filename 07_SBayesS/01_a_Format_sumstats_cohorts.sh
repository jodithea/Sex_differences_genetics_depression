#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=20GB
#PBS -l ncpus=1
#PBS -J 0-11

# This script formats input GWAS sumstats
# need to be in the GCTA-COJO .ma format = header row of SNP A1 A2 freq b se p N (SNP identifier (rsID), the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size)


### Environment ###




### Preamble ###

directory=/path/12_SBayesS_h2/

input_cohort=(/path/01_Format_sumstats/Females/*_MDD_female_sumstats_formatted_formetaanalysis.txt /path/01_Format_sumstats/Males/*_MDD_male_sumstats_formatted_formetaanalysis.txt)

name_cohort=$(echo ${input_cohort[${PBS_ARRAY_INDEX}]} | sed 's|/path/01_Format_sumstats/.*ales/||' | sed 's|_formatted_formetaanalysis.txt||')


### Run script ###

cd ${directory}

# Format cohort sumstats - all pre-formatted before meta-analysis to have same columns

awk 'NR==1 {print "SNP A1 A2 freq b se p N"} \
NR>1 {print $14, $3, $4, $5, $6, $7, $8, $9}' \
${input_cohort[${PBS_ARRAY_INDEX}]} \
> ${name_cohort}_formatted_SBayes.ma



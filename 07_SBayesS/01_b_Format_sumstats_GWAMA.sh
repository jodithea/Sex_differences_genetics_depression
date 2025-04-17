#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=20GB
#PBS -l ncpus=1
#PBS -J 0-1

# This script formats input GWAS sumstats
# need to be in the GCTA-COJO .ma format = header row of SNP A1 A2 freq b se p N (SNP identifier (rsID), the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size)


### Environment ###




### Preamble ###

directory=/path/12_SBayesS_h2/

input_GWAMA=(/path/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt /path/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt)

name_GWAMA=$(echo ${input_GWAMA[${PBS_ARRAY_INDEX}]} | sed 's|/path/03_Metal/.*ales/||' | sed 's|_QCed_rsID.txt||')



### Run script ###

cd ${directory}


# Format sex-stratified GWAS Meta-analysis sumstats

awk 'NR==1 {print "SNP A1 A2 freq b se p N"} \
NR>1 {print $20, toupper($2), toupper($3), $4, $8, $9, $10, $16}' \
${input_GWAMA[${PBS_ARRAY_INDEX}]} \
> ${name_GWAMA}_formatted_SBayes.ma


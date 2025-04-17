#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l mem=40GB
#PBS -l ncpus=1

# This script creates a Miami plot from two sets of GWAS summary statistics
# and adds rectangles highlighting evidence for segments detected in gwas-pw for model 1, 2, 3 or 4

### Environment ###

module load R/4.3.1

### Preamble ###

directory=/path/06_GWAS_pw/ # This is where the plots will be output

inputone=/path/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt # 1st set of summary stats - will be on top of Miami plot

inputtwo=/path/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt # 2nd set of summary stats - will be on bottom of Miami plot

nameone=Female # name to put as title for top of Miami plot

nametwo=Male # name to put as title for bottom of Miami plot

gwaspw_model1=GWAS_pw_femaleMDD_maleMDD_model1_ppa05.txt # file containing segments with evidence for model 1(association with first trait only)

gwaspw_model2=GWAS_pw_femaleMDD_maleMDD_model2_ppa05.txt # file containing segments with evidence for model 2 (association with 2nd trait only)

gwaspw_model3=GWAS_pw_femaleMDD_maleMDD_model3_ppa05.txt # file containing segments with evidence for model 3 (shared association with both traits)

gwaspw_model4=GWAS_pw_femaleMDD_maleMDD_model4_ppa05.txt # file containing segments with evidence for model 4 (two distinct associations with each trait)

p_threshold=5e-08 # p-value threshold for GWAS (where horizontal line will be placedon manhattan plot): 5e-08 for only common variants or 5e-09 for rare and common variants

p_threshold2=1e-06 # suggestive p-value threshold for GWAS (where 2nd lower horizontal line will be placed on manhattan plot): 1e-06 suggestive of significance

output=ppa05_femaleMDD_maleMDD # name of output (without file ending as .png will be added)

### Run script ###

cd ${directory}

Rscript --vanilla 06_Plot_Miami_gwaspw.R  ${inputone} ${inputtwo} ${nameone} ${nametwo} \
 ${gwaspw_model1} ${gwaspw_model2} ${gwaspw_model3} ${gwaspw_model4} \
 ${p_threshold} ${p_threshold2} \
 ${output}


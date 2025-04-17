#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l mem=10GB
#PBS -l ncpus=1

# This script pulls out which segments have evidence supporting each of the 4 models from gwas-pw
# 	model 1 [association only to phenotype 1 (maleMDD)] versus the null
#	model 2 [association only to phenotype 2 (maleBMI)] versus the null
#	model 3 [shared association to both phenotypes] versus the null
#	model 4 [two distinct associations, one to each phenotype] versus the null



### Environment ###

module load gzip/1.13


### Preamble ###

directory=/path/06_GWAS_pw/MDD_BMI/maleMDD_maleBMI_update/

input=GWAS_pw_maleMDD_maleBMI.segbfs.gz



### Submit script ###

cd ${directory}

gzip -dk ${input}


awk '($16 > 0.5)' ${input%.gz} > ${input%.segbfs.gz}_model1_ppa05.txt

awk '($17 > 0.5)' ${input%.gz} > ${input%.segbfs.gz}_model2_ppa05.txt

awk '($18 > 0.5)' ${input%.gz} > ${input%.segbfs.gz}_model3_ppa05.txt

awk '($19 > 0.5)' ${input%.gz} > ${input%.segbfs.gz}_model4_ppa05.txt

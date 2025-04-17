#!/bin/bash
#PBS -l walltime=00:20:00
#PBS -l mem=5GB
#PBS -l ncpus=1

# This script compares the shared regions found from gwas-pw 
# i.e. which regions are found in femaleMDD-femaleBMI model 3 AND maleMDD-maleBMI model3, only in femaleMDD-femaleBMI model 3, only in maleMDD-maleBMI model3


### Environment ###


### Preamble ###

directory=/path/06_GWAS_pw/MDD_BMI/

input1=/path/06_GWAS_pw/MDD_BMI/femaleMDD_femaleBMI/GWAS_pw_femaleMDD_femaleBMI_model3_ppa05.txt

input2=/path/06_GWAS_pw/MDD_BMI/maleMDD_maleBMI/GWAS_pw_maleMDD_maleBMI_model3_ppa05.txt

output=femaleMDD_femaleBMI_model3_vs_maleMDD_maleBMI_model3



### Submit script ###

cd ${directory}


header=$(head -n 1 ${input1})



echo ${header} > ${output}_regions_both.txt
awk 'NR==1 {next} FNR==1 {next} NR==FNR {seen[$1]; next} $1 in seen' ${input1} ${input2} >> ${output}_regions_both.txt


echo ${header} > ${output}_regions_female_only.txt
awk 'NR==1 {next} FNR==1 {next} NR==FNR {seen[$1]; next} !($1 in seen)' ${input2} ${input1} >> ${output}_regions_female_only.txt


echo ${header} > ${output}_regions_male_only.txt
awk 'NR==1 {next} FNR==1 {next} NR==FNR {seen[$1]; next} !($1 in seen)' ${input1} ${input2} >> ${output}_regions_male_only.txt




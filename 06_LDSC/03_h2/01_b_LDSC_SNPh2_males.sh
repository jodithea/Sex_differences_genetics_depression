#!/bin/bash
#PBS -l walltime=00:08:00
#PBS -l mem=5GB
#PBS -l ncpus=1

# This script uses LD score regression to check the intercept and calculate SNP-based heritability on the liability scale


### Environment ###

module load ldsc/20190815

### Preamble ###

# Edit below variables for your specific LDSC

directory=/path/04_LDSC/SNPh2/Males/

reference=/path/LDSC_reference/ # location of reference directory

file=/path/04_LDSC/Munge/Males/Metaanalysis_MDD_male_AllCohorts_LDSC_munged.sumstats.gz

name=Metaanalysis_MDD_male_AllCohorts

pop_prev=0.10 # lifetime prevalence of 10% in males

samp_prev=0.3290 # prevalence of trait in your sample (no. cases / (no. cases + no. controls))

# cases / (cases + controls) = 65,298/(65,298 + 132,724) = 0.3298



### Submit script ###

cd ${directory}

ldsc.py \
 --h2 ${file} \
 --ref-ld-chr ${reference}eur_w_ld_chr/ \
 --w-ld-chr ${reference}eur_w_ld_chr/ \
 --pop-prev ${pop_prev} \
 --samp-prev ${samp_prev} \
 --out ${name}_LDSC_h2


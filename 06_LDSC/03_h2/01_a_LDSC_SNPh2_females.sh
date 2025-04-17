#!/bin/bash
#PBS -l walltime=00:08:00
#PBS -l mem=5GB
#PBS -l ncpus=1

# This script uses LD score regression to check the intercept and calculate SNP-based heritability on the liability scale


### Environment ###

module load ldsc/20190815

### Preamble ###

# Edit below variables for your specific LDSC

directory=/path/04_LDSC/SNPh2/Females/

reference=/path/LDSC_reference/ # location of reference directory

file=/path/04_LDSC/Munge/Females/Metaanalysis_MDD_female_AllCohorts_LDSC_munged.sumstats.gz

name=Metaanalysis_MDD_female_AllCohorts

pop_prev=0.20 # lifetime prevalence of 20% in females

samp_prev=0.4499 # prevalence of trait in your sample (no. cases / (no. cases + no. controls))

# cases / (cases + controls) = 130,471 / (130,471 + 159,521) = 0.4499



### Submit script ###

cd ${directory}

ldsc.py \
 --h2 ${file} \
 --ref-ld-chr ${reference}eur_w_ld_chr/ \
 --w-ld-chr ${reference}eur_w_ld_chr/ \
 --pop-prev ${pop_prev} \
 --samp-prev ${samp_prev} \
 --out ${name}_LDSC_h2


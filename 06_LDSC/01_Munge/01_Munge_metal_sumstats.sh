#!/bin/bash
#PBS -l walltime=00:08:00
#PBS -l mem=5GB
#PBS -l ncpus=1
#PBS -J 0-6

# This script uses LDSC to munge the METAL meta-analysis results (all cohorts, and each of LOO)

### Environment ###

module load ldsc/20190815

### Preamble ###

# Edit below variables for your specific LDSC

directory=/path/04_LDSC/Munge/Females/

reference=/path/LDSC_reference/

files=(/path/03_Metal/Females/Metaanalysis_MDD_female_*_QCed_rsID.txt)

name=$(echo ${files[${PBS_ARRAY_INDEX}]} | sed 's|/path/03_Metal/Females/||' | sed 's|_QCed_rsID.txt||')


### Submit script ###

cd ${directory}

munge_sumstats.py \
 --sumstats ${files[${PBS_ARRAY_INDEX}]} \
 --merge-alleles ${reference}w_hm3.noMHC.snplist \
 --snp rsID_build37 \
 --ignore MarkerName \
 --N-col N_TOTAL \
 --a1 Allele1 \
 --a2 Allele2 \
 --p P \
 --signed-sumstats Effect,0 \
 --out ${name}_LDSC_munged

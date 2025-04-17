#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l mem=350GB
#PBS -l ncpus=15
#PBS -J 0-1

# This script uses SBayesS in GCTB to calculate SNP-based heritability, polygenicity and the selection parameter


### Environment ###

module load gctb/2.5.2


### Preamble ###

directory=/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/

LD_matrix=/reference/genepi/public_reference_panels/ldm_ukb_50k_bigset_2.8M/ukb50k_2.8M_shrunk_sparse.mldmlist #chr 1 - 22 doesnt include X chr

GWAS_sumstats=(${directory}*_formatted_SBayes.ma) # need to be in the GCTA-COJO .ma format = header row of SNP A1 A2 freq b se p N (SNP identifier (rsID), the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size)

name=$(basename "${GWAS_sumstats[${PBS_ARRAY_INDEX}]}" | sed 's|_formatted_SBayes.ma||')


### Submit script ###

cd ${directory}


gctb --sbayes S \
 --mldm ${LD_matrix} \
 --gwas-summary ${GWAS_sumstats[${PBS_ARRAY_INDEX}]} \
 --thread 15 \
 --chain-length 25000 \
 --burn-in 5000 \
 --num-chains 4 \
 --thin 10 \
 --seed 37904 \
 --write-mcmc-txt \
 --out ${name}_SBayesS

# *.parRes gives overall results
# Pi = polygenicity
# hsq = h2 = SNP-based heritability
# S = relationship between MAF and effect size (A negative relationship between effect size and MAF is a signature of negative (or purifying) selection)

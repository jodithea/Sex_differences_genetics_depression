#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l mem=15GB
#PBS -l ncpus=1

# This script creates Manhattan and QQ Plots 


### Environment ###

module load R/4.2.0


### Preamble ###

directory=/path/08_Sensitivity_analyses/Adams2024_topSNPs/ # This is where the plots will be output

p_threshold=7.17e-05 # p-value threshold for GWAS 0.05/697 independent SNPs = 7.17e-05 (based on Adams reporting 697 independent genome-wide sig SNPs)

file=Metaanalysis_MDD_GxS_fulldc_AllCohorts_QCed_rsID_Adams_sig_SNPs.txt

name=Metaanalysis_MDD_GxS_fulldc_AllCohorts_QCed_rsID_Adams_sig_SNPs



### Run script ###

cd ${directory}

Rscript --vanilla 03_Manhattan_plot.R ${file} ${name} ${p_threshold}


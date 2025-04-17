#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l mem=15GB
#PBS -l ncpus=1
#PBS -J 0-5

# This script creates Manhattan and QQ Plots from the GWAS summary statistics (using all SNPs)

### Environment ###

module load R/4.2.0


### Preamble ###

directory=/path/01_Format_sumstats/Females/ # This is where the plots will be output

p_threshold=5e-08 # p-value threshold for GWAS (where horizontal line will be placedon manhattan plot): 5e-08 for only common variants or 5e-09 for rare and common variants

files=(/path/01_Format_sumstats/Females/*_MDD_female_sumstats_formatted_formetaanalysis.txt)

name=$(echo ${files[${PBS_ARRAY_INDEX}]} | sed 's|/path/01_Format_sumstats/Females/||' | sed 's|_MDD_female_sumstats_formatted_formetaanalysis.txt||')


### Run script ###

cd ${directory}


# First change chr X for 23 so plots OK
awk '{ gsub("X", "23", $1) ; print }' ${files[${PBS_ARRAY_INDEX}]} > ${name}_temp.txt


Rscript --vanilla 04_Plots_Man_QQ.R ${name}_temp.txt ${name}_MDD_female_sumstats_formatted_formetaanalysis ${p_threshold} 

rm ${name}_temp.txt

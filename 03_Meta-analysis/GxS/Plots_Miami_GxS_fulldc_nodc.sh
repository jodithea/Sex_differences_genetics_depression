#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l mem=30GB
#PBS -l ncpus=1
#PBS -J 0-6

# This script creates Miami plot from two sets of GWAS summary statistics


### Environment ###

module load R/4.2.0


### Preamble ###

directory=/path/03_Metal/GxS/ # This is where the plots will be output

p_threshold=5e-08 # p-value threshold for GWAS (where horizontal line will be placedon manhattan plot): 5e-08 for only common variants or 5e-09 for rare and common variants

p_threshold2=1e-06

inputone=(/path/03_Metal/GxS/fulldc/*_QCed_rsID.txt) # name of your first input file - will be on top of Miami plot, e.g. GWAS.fastGWA

inputtwo=(/path/03_Metal/GxS/nodc/*_QCed_rsID.txt) # name of your second input file - will be on bottom of Miami plot, e.g. GWAS2.fastGWA

labelone=Full_DC # name to put as title for top of Miami plot

labeltwo=No_DC # name to put as title for bottom of Miami plot

name=$(echo ${inputone[${PBS_ARRAY_INDEX}]} | sed 's|/path/03_Metal/GxS/fulldc/||' | sed 's|_fulldc||' | sed 's|_QCed_rsID.txt||') # name of output (without file ending as .png will be added)

### Run script ###

cd ${directory}


Rscript --vanilla Plots_Miami_2siglines.R \
 ${inputone[${PBS_ARRAY_INDEX}]} \
 ${inputtwo[${PBS_ARRAY_INDEX}]} \
 ${labelone} \
 ${labeltwo} \
 ${name}_fulldc_nodc_QCed_rsID \
 ${p_threshold} \
 ${p_threshold2}


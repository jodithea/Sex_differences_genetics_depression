#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l mem=5GB
#PBS -l ncpus=1

# This script tests whether rg of female MDD to male MDD is significantly different to 1
# And if the rg between female MDD and a trait is significantly different to the rg between male MDD and the same trait (using Z-score method which is not theoretically appropriate but commonly done)
# Used benjamini hochberg method for multiple comparisons

### Environment ###

module load R/4.3.1


### Run script ###

cd /path/04_LDSC/SNPrg/

Rscript --vanilla 01_Compare_female_male_rg_traits_Z_score_method.R

#!/bin/bash
#PBS -l walltime=00:08:00
#PBS -l mem=5GB
#PBS -l ncpus=1
#PBS -J 0-6

# This script uses LD score regression to check the intercept of the meta-analysis results


### Environment ###

module load ldsc/20190815

### Preamble ###

directory=/path/04_LDSC/Intercept/Females/

reference=/path/LDSC_reference/ # location of reference directory

files=(/path/04_LDSC/Munge/Females/*_LDSC_munged.sumstats.gz)

name=$(echo ${files[${PBS_ARRAY_INDEX}]} | sed 's|/path/04_LDSC/Munge/Females/||' | sed 's|_LDSC_munged.sumstats.gz||')


### Submit script ###

cd ${directory}

ldsc.py \
 --h2 ${files[${PBS_ARRAY_INDEX}]} \
 --ref-ld-chr ${reference}eur_w_ld_chr/ \
 --w-ld-chr ${reference}eur_w_ld_chr/ \
 --out ${name}_LDSC_h2



# Create output with info on intercept

echo "LDSC Intercept and SD:" >> ${name}_LDSCintercept.txt
grep "Intercept:" ${name}_LDSC_h2.log >> ${name}_LDSCintercept.txt

intercept=$(grep "Intercept:" ${name}_LDSC_h2.log | awk '{print $2}')
sd=$(grep "Intercept:" ${name}_LDSC_h2.log | awk '{print $3}' | sed 's/[()]//g')
lower_ci=$(echo "$intercept - (1.96 * $sd)" | bc)
upper_ci=$(echo "$intercept + (1.96 * $sd)" | bc)
echo "LDSC Intercept 95% CI: Lower = $lower_ci, Upper = $upper_ci" >> ${name}_LDSCintercept.txt

# In *_LDSC_h2.log file the intercept and the SD (in brackets) is reported
# Calculate the intercept's 95% CI
#  95% CI = 1.96 * SD
#  Therefore, intercept +/- (SD x 1.96) gives 95% CI
# Intercept not sig dif from 1 suggests lambda inflation is true polygenicity and not due to confounding.
# Intercept sig > 1 suggests confounding
# Intercept sig < 1 suggests overcorrection
# Also check the ratio - expect to be not significantly different from 0

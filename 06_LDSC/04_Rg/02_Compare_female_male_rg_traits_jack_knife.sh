#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l mem=10GB
#PBS -l ncpus=1
#PBS -J 0-17


# This script compares the rg of female MDD with a trait and rg of male MDD with the same trait using jack knife method


### Environment ###

module load R/4.3.1


### Preamble ###

directory=/path/04_LDSC/SNPrg/jack_knife/

infile_rg1=/path/04_LDSC/SNPrg/Females/rg_femaleMDD_correlation_table.txt

infile_rg2=/path/04_LDSC/SNPrg/Males/rg_maleMDD_correlation_table.txt


trait=(
  "AdamsMDD_sumstats"
  "BloklandMDDfemale_sumstats"
  "BloklandMDDmale_sumstats"
  "SilveiraMDDfemale_sumstats"
  "SilveiraMDDmale_sumstats"
  "ADHD_sumstats"
  "Anxiety_sumstats"
  "Bipolar_sumstats"
  "PTSD_sumstats"
  "Schizophrenia_sumstats"
  "BMI_sumstats"
  "Waist_to_hip_ratio_sumstats"
  "EA_sumstats"
  "Metabolic_syndrome_sumstats"
  "Drinks_per_week_sumstats"
  "Ever_smoked_regularly_sumstats"
  "BMI_female_sumstats"
  "BMI_male_sumstats"
)

# Corresponding identifiers to search in filenames 
trait_files=(
  "Adams_Full_EUR_2025"
  "meta_STDERR_mdd_eur_auto_xchr_F1"
  "meta_STDERR_mdd_eur_auto_xchr_M1"
  "GWAS.BroadMDD_Females"
  "GWAS.BroadMDD_Males"
  "ADHD2022_iPSYCH_deCODE_PGC"
  "ANX_EUR.txt.gz"
  "Bipolar_Neff"
  "PTSD_Neff"
  "Schizophrenia_Neff"
  "Meta-analysis_Locke_et_al_plus_UKBiobank"
  "fat-distn.giant.ukbb.meta-analysis.whradjbmi.combined"
  "EA4_additive_excl_23andMe"
  "GCST90444487"
  "GSCAN_DrnkWk_2022_GWAS_SUMMARY_STATS_EUR"
  "GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR"
  "bmi.giant-ukbb.meta-analysis.females.23May2018"
  "bmi.giant-ukbb.meta-analysis.males.23May2018"
)

# Build file arrays
female_dir=/path/04_LDSC/SNPrg/Females/

male_dir=/path/04_LDSC/SNPrg/Males/

for i in "${!trait[@]}"; do
  t="${trait_files[$i]}"
  infile_c1[$i]=$(ls ${female_dir}/*${t}*.gencov.delete)
  infile_c2[$i]=$(ls ${male_dir}/*${t}*.gencov.delete)

  infile_hsq1[$i]="${infile_c1[$i]/gencov/hsq1}"
  infile_hsq2[$i]="${infile_c1[$i]/gencov/hsq2}"
  infile_hsq3[$i]="${infile_c2[$i]/gencov/hsq1}"
  infile_hsq4[$i]="${infile_c2[$i]/gencov/hsq2}"
done



### Run script ###

cd ${directory}

echo "Running job index: ${PBS_ARRAY_INDEX}"
echo "Trait: ${trait[${PBS_ARRAY_INDEX}]}"
echo "Using infile_c1: ${infile_c1[${PBS_ARRAY_INDEX}]}"
echo "Using infile_c2: ${infile_c2[${PBS_ARRAY_INDEX}]}"
echo "Using infile_hsq1: ${infile_hsq1[${PBS_ARRAY_INDEX}]}"
echo "Using infile_hsq2: ${infile_hsq2[${PBS_ARRAY_INDEX}]}"
echo "Using infile_hsq3: ${infile_hsq3[${PBS_ARRAY_INDEX}]}"
echo "Using infile_hsq4: ${infile_hsq4[${PBS_ARRAY_INDEX}]}"

Rscript --vanilla 02_Compare_female_male_rg_traits_jack_knife.R \
	${infile_rg1} \
	${infile_rg2} \
	${infile_c1[${PBS_ARRAY_INDEX}]} \
	${infile_c2[${PBS_ARRAY_INDEX}]} \
	${infile_hsq1[${PBS_ARRAY_INDEX}]} \
	${infile_hsq2[${PBS_ARRAY_INDEX}]} \
	${infile_hsq3[${PBS_ARRAY_INDEX}]} \
	${infile_hsq4[${PBS_ARRAY_INDEX}]} \
	${trait[${PBS_ARRAY_INDEX}]} \
	${directory}

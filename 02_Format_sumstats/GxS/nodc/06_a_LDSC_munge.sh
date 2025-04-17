#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l mem=5GB
#PBS -l ncpus=1


# This script uses LD score regression to format summary stats with munge - for gene by environment results
# Then the munged summary stats go in to calculate the LDSC intercept

### Environment ###

module load ldsc/20190815

### Preamble ###

# Edit below variables for your specific LDSC

directory=/path/01_Format_sumstats/GxS/nodc/ # working directory

reference=/path/LDSC_reference/ # location of reference directory



### Submit script ###

cd ${directory}

# rather than doing in parallel just do each manually
# They don't take long to run and there are differences depending on if the input sumstats contain the imputation score
# From LDSC docs: Imputation quality is correlated with LD Score, and low imputation quality yields lower test statistics, so imputation quality is a confounder for LD Score regression. To prevent bias from variable imputation quality, we usually remove poorly-imputed SNPs by filtering on INFO > 0.9. The scz and bip summary statistics that we're using for this tutorial have INFO columns, so munge_sumstats.py will automatically perform the filtering. If you're using summary statistics that don't come with an INFO column, we recommend filtering to HapMap3 SNPs (using the --merge or --merge-alleles flags), because these seem to be well-imputed in most studies.
# For All Of Us set info to 1 for all SNPs as was whole genome sequenced
# Blokland21 has no info column so use advice to restrict to hapmap SNPs and don't  filter on info


mkdir Munged_sumstats

# AGDS

munge_sumstats.py \
 --sumstats AGDS_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed.txt \
 --merge-alleles ${reference}w_hm3.noMHC.snplist \
 --snp rsID_build37 \
 --frq FREQA1 \
 --signed-sumstats Z,0 \
 --p P_G_by_E \
 --ignore MARKER_build37,MAF \
 --out Munged_sumstats/AGDS_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged



# Bionic

munge_sumstats.py \
 --sumstats Bionic_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed.txt \
 --merge-alleles ${reference}w_hm3.noMHC.snplist \
 --snp rsID_build37 \
 --frq FREQA1 \
 --signed-sumstats Z,0 \
 --p P_G_by_E \
 --ignore MARKER_build37,MAF \
 --out Munged_sumstats/Bionic_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged


# UKB

munge_sumstats.py \
 --sumstats UKB_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed.txt \
 --merge-alleles ${reference}w_hm3.noMHC.snplist \
 --snp rsID_build37 \
 --frq FREQA1 \
 --signed-sumstats Z,0 \
 --p P_G_by_E \
 --ignore MARKER_build37,MAF \
 --out Munged_sumstats/UKB_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged


# GLAD

munge_sumstats.py \
 --sumstats GLAD_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed.txt \
 --merge-alleles ${reference}w_hm3.noMHC.snplist \
 --snp rsID_build37 \
 --frq FREQA1 \
 --signed-sumstats Z,0 \
 --p P_G_by_E \
 --ignore MARKER_build37,MAF \
 --out Munged_sumstats/GLAD_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged


# Blokland21

munge_sumstats.py \
 --sumstats Blokland21_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed.txt \
 --merge-alleles ${reference}w_hm3.noMHC.snplist \
 --snp rsID_build37 \
 --frq FREQA1 \
 --signed-sumstats BETA_G_by_E,0 \
 --p P_G_by_E \
 --ignore MARKER_build37,INFO,MAF \
 --out Munged_sumstats/Blokland21_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged


# All Of Us

awk '{ gsub("NA", "1", $15) ; print }' AllOfUs_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed.txt > temp_AllOfUs.txt

munge_sumstats.py \
 --sumstats temp_AllOfUs.txt \
 --merge-alleles ${reference}w_hm3.noMHC.snplist \
 --snp rsID_build37 \
 --frq FREQA1 \
 --signed-sumstats Z,0 \
 --p P_G_by_E \
 --ignore MARKER_build37,MAF \
 --out Munged_sumstats/AllOfUs_MDD_GxS_nodc_sumstats_formatted_formetaanalysis_transformed_LDSC_munged

rm temp_AllOfUs.txt

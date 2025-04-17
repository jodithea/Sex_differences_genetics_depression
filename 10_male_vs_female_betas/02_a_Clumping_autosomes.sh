#!/bin/bash
#PBS -l walltime=4:00:00
#PBS -l mem=20GB
#PBS -l ncpus=5
#PBS -J 1-22

# This script used plink to carry out clumping to identify the number of independent SNPs that are genome-wide significant

### Preamble ###

directory=/path/09_Effect_sizes_plots/Adams_SNPs/

chr=${PBS_ARRAY_INDEX}

input=Adams_Full_EUR_2025_FromIMB_formatted.txt


### Environment ###
module load plink/1.90b6.8


### Submit script ###
cd ${directory}

# Autosomes, using rsID to match bfile to sumstats

plink --bfile /path/Meta-analysis/1000G_Phase3_v5/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_rsID_no_missing_dups \
  --keep /path/Meta-analysis/1000G_Phase3_v5/EUR_population_IDs.txt \
  --clump ${input} \
  --clump-p1 0.00000005 \
  --clump-p2 1 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --out ${directory}/clumped_Adams_SNPs_chr${chr}

# Uses LD ref panel = 1000 Genomes Phase 3 files with marker name of the format rsID according to the marker names in NCBI dbSNP build 155 human geome build 37 patch 13

# -- clump                     text files with a header line, a column containing variant IDs (SNP), and another column containing p-values (P)
# --clump-p1 0.00000005        Significance threshold for index SNPs (i.e. the most sig SNP in each clump must have a p-value at this value or lower) = 5e-08 so only create clumps for genome-wide sig SNPs
# --clump-p2 1                 Secondary significance threshold for clumped SNPs (i.e. SNPs with p-value at this value or lower will be listed as a SNP within the appropriate clump - using 1 means all SNPs are assigned into a clump)
# --clump-r2 0.1               LD threshold for clumping based on maximum likelihood haplotype frequency estimates
# --clump-kb 1000              Physical distance threshold for clumping (The maximum distance from the lead variant for SNPs to be considered to be clumped with it)



# OUPUT
#      CHR     Chromosome code
#      F       Results fileset code (1,2,...)
#      SNP     SNP identifier
#      BP      Physical position of SNP (base-pairs)
#      TOTAL   Total number of other SNPs in clump (i.e. passing --clump-kb and --clump-r2 thresholds)
#      NSIG    Number of clumped SNPs that are not significant ( p > 0.05 )
#      S05     Number of clumped SNPs 0.01 < p < 0.05
#      S01     Number of clumped SNPs 0.001 < p < 0.01
#      S001    Number of clumped SNPs 0.0001 < p < 0.001
#      S0001   Number of clumped SNPs p < 0.0001
#      SP2     List of SNPs names (and fileset code) clumped and significant at --clump-p2

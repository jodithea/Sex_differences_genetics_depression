#!/bin/bash
#PBS -l walltime=2:00:00
#PBS -l mem=10GB
#PBS -l ncpus=5

# This script used plink to carry out clumping to identify the number of independent SNPs that are genome-wide significant

### Preamble ###

directory=/path/05_Clumping/Females/

input=Metaanalysis_MDD_female_AllCohorts_QCed_rsID_SNPs.txt

### Environment ###
module load plink/1.90b6.8


### Submit script ###
cd ${directory}


# X chromosome
plink --bfile /path/Meta-analysis/1000G_Phase3_v5/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes_CHRBPA1A2_no_missing_dups \
  --keep /path/Meta-analysis/1000G_Phase3_v5/EUR_population_IDs.txt \
  --clump ${input} \
  --clump-p1 0.0001 \
  --clump-p2 1 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --out ${directory}/clumped_female_chrX

# Uses LD ref panel = 1000 Genomes Phase 3 files with marker name of the format CHR:BP:A1:A2 according to the marker names in NCBI dbSNP build 155 human geome build 37 patch 13
# Before the meta-analysis I formatted the SNPs to use a consistent marker name using the format CHR:BP:A1:A2 according to the marker names in NCBI dbSNP build 155 human geome build 37 patch 13
# So each SNP in the meta-analysis input file will match the LD ref panel

# -- clump                     text files with a header line, a column containing variant IDs (SNP), and another column containing p-values (P)
# --clump-p1 0.0001            Significance threshold for index SNPs (i.e. the most sig SNP in each clump must have a p-value at this value or lower)
# --clump-p2 1                 Secondary significance threshold for clumped SNPs (i.e. SNPs with p-value at this value or lower will be listed as a SNP within the appropriate clump - using 1 means all SNPs are assigned into a clump)
# --clump-r2 0.50              LD threshold for clumping based on maximum likelihood haplotype frequency estimates
# --clump-kb 250               Physical distance threshold for clumping (The maximum distance from the lead variant for SNPs to be considered to be clumped with it)



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

directory="/path/06_GWAS_pw/MDD_BMI/femaleMDD_femaleBMI_update/" # working directory

femaleMDD="/path/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt"

femaleBMI="/path/Summary_statistics_published/BMI_sex_Pulit_2018/bmi.giant-ukbb.meta-analysis.females.23May2018.txt.gz"


### LOAD PACKAGES ###
library(tidyverse)

### READ IN DATA ###
femaleMDD_df <- read.table(femaleMDD, header = T, stringsAsFactors = F)

femaleBMI_df <- read.table(gzfile(femaleBMI), header = T, stringsAsFactors = F)


### FORMAT FOR FGWAS AND GWAS-PW ###

## Create Z-score and Variance columns ##

# Z = the signed Z score measuring the evidence for association to phenotype at the SNP = Effect/StdErr
# V = variance = SE^2
# Add "chr" before chromosome number (so matches style in bed file input for fgwas and gwas-pw)
# MDD - change A1 and A2 to upper case
# MDD MarkerName based on refernece dbSNP file to create consistent for meta-analysis but not actually reflecting A1 then A2 (markername could flip alleles). So make markername that can be used for matching to ensure same alleles

femaleMDD_format_df <- femaleMDD_df %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2),
         MarkerName_match = paste(CHR, BP, Allele1, Allele2, sep = ":"),
         Z = Effect / StdErr,
         V = StdErr * StdErr,
         CHR = paste("chr", CHR, sep = "")) %>%
  rename(SE = StdErr,
         rsID = rsID_build37,
         FREQA1 = Freq1,
         N = N_TOTAL) %>%
  select(MarkerName_match, rsID, CHR, BP, Allele1, Allele2, Effect, Z, V, SE, P, FREQA1, N) %>%
  rename_with(~ paste0(., "_femaleMDD"))


femaleBMI_format_df <- femaleBMI_df %>%
  mutate(Z = BETA / SE,
         V = SE * SE,
         MarkerName = paste(CHR, POS, Tested_Allele, Other_Allele, sep = ":"),
         MarkerNameflipped = paste(CHR, POS, Other_Allele, Tested_Allele, sep = ":"),
         CHR = paste("chr", CHR, sep = ""),
         rsID = str_extract(SNP, "[^:]+")) %>%
  rename(Allele1 = Tested_Allele,
         Allele2 = Other_Allele,
         BP = POS,
         FREQA1 = Freq_Tested_Allele,
         Effect = BETA) %>%
  select(MarkerName, MarkerNameflipped, rsID, CHR, BP, Allele1, Allele2, Effect, Z, V, SE, P, FREQA1, N) %>%
  rename_with(~ paste0(., "_femaleBMI"))



## Compare SNPs ##

# fgwas matches to bed file that only includes range of positions (no alleles)
# so need to make sure that when SNPs are matched across sumstats that the SNPs are the same
# First check that the tested allele is the same in each sumstats
# If tested allele is swapped (i.e. a1 MDD = a2 BMI) then need to flip sign of effect size so effect size can be compared across sumstats (as Z = effect / SE this would mean flipping sign of Z)

merged_df <- femaleMDD_format_df %>%
  inner_join(femaleBMI_format_df, by = c("MarkerName_match_femaleMDD" = "MarkerName_femaleBMI")) %>%
  rename(MarkerName = MarkerName_match_femaleMDD)

merged_flipped_df <- femaleMDD_format_df %>%
  inner_join(femaleBMI_format_df, by = c("MarkerName_match_femaleMDD" = "MarkerNameflipped_femaleBMI")) %>%
  rename(MarkerName = MarkerName_match_femaleMDD)

merged_flipped_df <- merged_flipped_df %>%
  mutate(Z_femaleBMI = -1 * Z_femaleBMI)

femaleMDD_femaleBMI_df <- merged_df %>%
  bind_rows(merged_flipped_df)

# Sort rows by chromosomal position
femaleMDD_femaleBMI_df <- femaleMDD_femaleBMI_df %>%
  arrange(CHR_femaleMDD, BP_femaleMDD)

# Remove chr 23 as not in bed file
femaleMDD_femaleBMI_df <- femaleMDD_femaleBMI_df %>%
  filter(CHR_femaleMDD != 23)

# fgwas can only deal with one SNP per position (can't deal with multiallelic SNPs)
# So need to only keep one SNP per positions - make sure this is the same in both sumstats so not comparing different SNPs at the same position across the two sumstats
# in combined file, when there are multiple SNPs with same MarkerName (i.e. multiallelic site) keep first SNP

femaleMDD_femaleBMI_distinct_df <- femaleMDD_femaleBMI_df %>%
  distinct(CHR_femaleMDD, BP_femaleMDD, .keep_all = T)

# As a santiy check - view all multiallelic sites
# duplicates_df <- femaleMDD_femaleBMI_df %>%
#   group_by(CHR_femaleMDD, BP_femaleMDD) %>%
#   filter(n() > 1) %>%
#   ungroup()

## fgwas files ##

# Split combined file and create one for each sumstats as input for fgwas: SNPID, CHR, POS, Z, SE, F, N
# Documentation states if include SE 'this column will be taken as a direct estimate of the standard errors of the regression coefficient and will override the sample size and allele frequency columns.'
# However, if include SE and no F or N then get WARNING: detected SE in header, will override F and N and ERROR: cannot find F in header
# So include F = : the allele frequency of one of the alleles of the SNP and N (for case/control studies should provide NCASE and NCONTROL but SE should override this anyway)
femaleMDD_formatted_fgwas <- femaleMDD_femaleBMI_distinct_df %>%
  select(MarkerName, CHR_femaleMDD, BP_femaleMDD, Z_femaleMDD, SE_femaleMDD, FREQA1_femaleMDD, N_femaleMDD) %>%
  rename(SNPID = MarkerName,
         CHR = CHR_femaleMDD,
         POS = BP_femaleMDD,
         Z = Z_femaleMDD,
         SE = SE_femaleMDD,
         F = FREQA1_femaleMDD,
         N = N_femaleMDD)

outfile <- paste(directory, "female_MDD_sumstats_formatted_fgwas.txt", sep="")
write.table(femaleMDD_formatted_fgwas, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


femaleBMI_formatted_fgwas <- femaleMDD_femaleBMI_distinct_df %>%
  select(MarkerName, CHR_femaleBMI, BP_femaleBMI, Z_femaleBMI, SE_femaleBMI, FREQA1_femaleBMI, N_femaleBMI) %>%
  rename(SNPID = MarkerName,
         CHR = CHR_femaleBMI,
         POS = BP_femaleBMI,
         Z = Z_femaleBMI,
         SE = SE_femaleBMI,
         F = FREQA1_femaleBMI,
         N = N_femaleBMI)

outfile <- paste(directory, "female_BMI_sumstats_formatted_fgwas.txt", sep="")
write.table(femaleBMI_formatted_fgwas, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


## gwas-pw file ##

# Format combined file so can go into gwas-pw: SNPID, CHR, POS, Z[sumstats1], V[sumstats1], Z[sumstats2], V[sumstats2]
femaleMDD_femaleBMI_formatted_gwaspw <- femaleMDD_femaleBMI_distinct_df %>%
  select(MarkerName, CHR_femaleMDD, BP_femaleMDD, Z_femaleMDD, V_femaleMDD, Z_femaleBMI, V_femaleBMI) %>%
  rename(SNPID = MarkerName,
         CHR = CHR_femaleMDD,
         POS = BP_femaleMDD)

outfile <- paste(directory, "femaleMDD_femaleBMI_sumstats_formatted_gwaspw.txt", sep="")
write.table(femaleMDD_femaleBMI_formatted_gwaspw, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



## Also save formatted sumstats ##
# Use femaleMDD_femaleBMI_distinct_df (two sumstats combined and SNPs matched appropriately if alleles flipped, top SNP kept for multialleleic sites, ordered by chr and chr23 (X) removed)
# Use to match fgwas results to this file when calculating correlation

femaleMDD_sumstats <- femaleMDD_femaleBMI_distinct_df %>%
  select(MarkerName, ends_with("femaleMDD")) %>%
  rename_with(~ str_replace_all(., "_femaleMDD", ""))

outfile <- paste(directory, "femaleMDD_sumstats_matched_distinct.txt", sep="")
write.table(femaleMDD_sumstats, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


femaleBMI_sumstats <- femaleMDD_femaleBMI_distinct_df %>%
  select(MarkerName, ends_with("femaleBMI"), -MarkerNameflipped_femaleBMI, -MarkerName_femaleBMI) %>%
  rename_with(~ str_replace_all(., "_femaleBMI", ""))

outfile <- paste(directory, "femaleBMI_sumstats_matched_distinct.txt", sep="")
write.table(femaleBMI_sumstats, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


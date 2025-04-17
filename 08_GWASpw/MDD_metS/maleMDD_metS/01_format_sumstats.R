directory="/path/06_GWAS_pw/MDD_metS/maleMDD_metS/" # working directory

maleMDD="/path/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt"

metS="/path/Summary_statistics_published/Metabolic_syndrome_Park_2024/GCST90444487.tsv"



### LOAD PACKAGES ###
library(tidyverse)

### READ IN DATA ###
maleMDD_df <- read.table(maleMDD, header = T, stringsAsFactors = F)

metS_df <- read.table(metS, header = T, stringsAsFactors = F)


### FORMAT FOR FGWAS AND GWAS-PW ###

## Create Z-score and Variance columns ##

# Z = the signed Z score measuring the evidence for association to phenotype at the SNP = Effect/StdErr
# V = variance = SE^2
# Add "chr" before chromosome number (so matches style in bed file input for fgwas and gwas-pw)
# MDD - change A1 and A2 to upper case
# MDD MarkerName based on refernece dbSNP file to create consistent for meta-analysis but not actually reflecting A1 then A2 (markername could flip alleles). So make markername that can be used for matching to ensure same alleles

maleMDD_format_df <- maleMDD_df %>%
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
  rename_with(~ paste0(., "_maleMDD"))


metS_format_df <- metS_df %>%
  mutate(Z = beta / standard_error,
         V = standard_error * standard_error,
         MarkerName = paste(chromosome, base_pair_location, effect_allele, other_allele, sep = ":"),
         MarkerNameflipped = paste(chromosome, base_pair_location, other_allele, effect_allele, sep = ":"),
         CHR = paste("chr", chromosome, sep = ""),
         rsID = str_extract(rs_id, "[^:]+"),
         N = 1384348) %>%
  rename(Allele1 = effect_allele,
         Allele2 = other_allele,
         BP = base_pair_location,
         FREQA1 = effect_allele_frequency,
         Effect = beta,
         SE = standard_error,
         P = p_value) %>%
  select(MarkerName, MarkerNameflipped, rsID, CHR, BP, Allele1, Allele2, Effect, Z, V, SE, P, FREQA1, N) %>%
  rename_with(~ paste0(., "_metS"))



## Compare SNPs ##

# fgwas matches to bed file that only includes range of positions (no alleles)
# so need to make sure that when SNPs are matched across sumstats that the SNPs are the same
# First check that the tested allele is the same in each sumstats
# If tested allele is swapped (i.e. a1 MDD = a2 BMI) then need to flip sign of effect size so effect size can be compared across sumstats (as Z = effect / SE this would mean flipping sign of Z)

merged_df <- maleMDD_format_df %>%
  inner_join(metS_format_df, by = c("MarkerName_match_maleMDD" = "MarkerName_metS")) %>%
  rename(MarkerName = MarkerName_match_maleMDD)

merged_flipped_df <- maleMDD_format_df %>%
  inner_join(metS_format_df, by = c("MarkerName_match_maleMDD" = "MarkerNameflipped_metS")) %>%
  rename(MarkerName = MarkerName_match_maleMDD)

merged_flipped_df <- merged_flipped_df %>%
  mutate(Z_metS = -1 * Z_metS)

maleMDD_metS_df <- merged_df %>%
  bind_rows(merged_flipped_df)

# Sort rows by chromosomal position
maleMDD_metS_df <- maleMDD_metS_df %>%
  arrange(CHR_maleMDD, BP_maleMDD)

# Remove chr 23 as not in bed file
maleMDD_metS_df <- maleMDD_metS_df %>%
  filter(CHR_maleMDD != 23)

# fgwas can only deal with one SNP per position (can't deal with multiallelic SNPs)
# So need to only keep one SNP per positions - make sure this is the same in both sumstats so not comparing different SNPs at the same position across the two sumstats
# in combined file, when there are multiple SNPs with same MarkerName (i.e. multiallelic site) keep first SNP

maleMDD_metS_distinct_df <- maleMDD_metS_df %>%
  distinct(CHR_maleMDD, BP_maleMDD, .keep_all = T)

# As a santiy check - view all multiallelic sites
# duplicates_df <- maleMDD_metS_df %>%
#   group_by(CHR_maleMDD, BP_maleMDD) %>%
#   filter(n() > 1) %>%
#   ungroup()

## fgwas files ##

# Split combined file and create one for each sumstats as input for fgwas: SNPID, CHR, POS, Z, SE, F, N
# Documentation states if include SE 'this column will be taken as a direct estimate of the standard errors of the regression coefficient and will override the sample size and allele frequency columns.'
# However, if include SE and no F or N then get WARNING: detected SE in header, will override F and N and ERROR: cannot find F in header
# So include F = : the allele frequency of one of the alleles of the SNP and N (for case/control studies should provide NCASE and NCONTROL but SE should override this anyway)
maleMDD_formatted_fgwas <- maleMDD_metS_distinct_df %>%
  select(MarkerName, CHR_maleMDD, BP_maleMDD, Z_maleMDD, SE_maleMDD, FREQA1_maleMDD, N_maleMDD) %>%
  rename(SNPID = MarkerName,
         CHR = CHR_maleMDD,
         POS = BP_maleMDD,
         Z = Z_maleMDD,
         SE = SE_maleMDD,
         F = FREQA1_maleMDD,
         N = N_maleMDD)

outfile <- paste(directory, "male_MDD_sumstats_formatted_fgwas.txt", sep="")
write.table(maleMDD_formatted_fgwas, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


metS_formatted_fgwas <- maleMDD_metS_distinct_df %>%
  select(MarkerName, CHR_metS, BP_metS, Z_metS, SE_metS, FREQA1_metS, N_metS) %>%
  rename(SNPID = MarkerName,
         CHR = CHR_metS,
         POS = BP_metS,
         Z = Z_metS,
         SE = SE_metS,
         F = FREQA1_metS,
         N = N_metS)

outfile <- paste(directory, "metS_sumstats_formatted_fgwas.txt", sep="")
write.table(metS_formatted_fgwas, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


## gwas-pw file ##

# Format combined file so can go into gwas-pw: SNPID, CHR, POS, Z[sumstats1], V[sumstats1], Z[sumstats2], V[sumstats2]
maleMDD_metS_formatted_gwaspw <- maleMDD_metS_distinct_df %>%
  select(MarkerName, CHR_maleMDD, BP_maleMDD, Z_maleMDD, V_maleMDD, Z_metS, V_metS) %>%
  rename(SNPID = MarkerName,
         CHR = CHR_maleMDD,
         POS = BP_maleMDD)

outfile <- paste(directory, "maleMDD_metS_sumstats_formatted_gwaspw.txt", sep="")
write.table(maleMDD_metS_formatted_gwaspw, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



## Also save formatted sumstats ##
# Use maleMDD_metS_distinct_df (two sumstats combined and SNPs matched appropriately if alleles flipped, top SNP kept for multialleleic sites, ordered by chr and chr23 (X) removed)
# Use to match fgwas results to this file when calculating correlation

maleMDD_sumstats <- maleMDD_metS_distinct_df %>%
  select(MarkerName, ends_with("maleMDD")) %>%
  rename_with(~ str_replace_all(., "_maleMDD", ""))

outfile <- paste(directory, "maleMDD_sumstats_matched_distinct.txt", sep="")
write.table(maleMDD_sumstats, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


metS_sumstats <- maleMDD_metS_distinct_df %>%
  select(MarkerName, ends_with("metS"), -MarkerNameflipped_metS, -MarkerName_metS) %>%
  rename_with(~ str_replace_all(., "_metS", ""))

outfile <- paste(directory, "metS_sumstats_matched_distinct.txt", sep="")
write.table(metS_sumstats, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



# This script compares genome-wide sig SNPs from sex-stratified meta-analysis to the same SNPs in our replication GWAS done in Generation Scotland (with individuals iverlapping with UK Biobank removed)

directory = "/path/11_Replication_GWAS/"

femaleGWAS_MA = "/path/Meta-analysis/Freeze2/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt"

maleGWAS_MA = "/path/Meta-analysis/Freeze2/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt"

clumpfile_female = "/path/05_Clumping/Females/clumped_female_all_results.txt"

clumpfile_male = "/path/05_Clumping/Males/clumped_male_all_results.txt"

femaleGWAS_GenScot = "/path/GenerationScotland_Poppy_results/GenerationScotland_SexMDDGWAS/GS_MDD_females_allchrs.fastGWA"

maleGWAS_GenScot = "/path/GenerationScotland_Poppy_results/GenerationScotland_SexMDDGWAS/GS_MDD_males_allchrs.fastGWA"

#########################################################################

### LOAD PACKAGES ###

library(tidyverse)



### READ IN DATA ###

femaleGWAS_MA_df <- read.table(femaleGWAS_MA, header = TRUE, stringsAsFactors = F)

maleGWAS_MA_df <- read.table(maleGWAS_MA, header = TRUE, stringsAsFactors = F)

clump_df_female <- read.table(clumpfile_female, header = TRUE, stringsAsFactors = F)

clump_df_male <- read.table(clumpfile_male, header = TRUE, stringsAsFactors = F)

femaleGWAS_GenScot_df <- read.table(femaleGWAS_GenScot, header = TRUE, stringsAsFactors = F)

maleGWAS_GenScot_df <- read.table(maleGWAS_GenScot, header = TRUE, stringsAsFactors = F)


#########FEMALES##############

### EXTRACT GENOME-WIDE SIG SNPs ###

# Extract only lead independent sig SNPs and match with GenScot sumstats (ensuring alleles matched correctly and if flipped then flip beta)

female_ind_sig_SNPs <- clump_df_female %>%
  filter(P < 5e-08) %>%
  select(SNP) %>%
  pull(SNP)

# Pull out lead ind sig SNPs using MarkerName (same MarkerName used in sumstats and clumping)
# MarkerName is based on dbSNP ref and doesn't necessarily match A1:A2 so make new markername to match with GenScot replication sumstats
femaleSNPs_MA_df <- femaleGWAS_MA_df %>%
  filter(MarkerName %in% female_ind_sig_SNPs) %>%
  rename(SE = StdErr,
         rsID = rsID_build37,
         FREQA1 = Freq1) %>%
  mutate(MarkerName_match = paste(CHR, BP, toupper(Allele1), toupper(Allele2), sep = ":")) %>%
  select(MarkerName_match, rsID, CHR, BP, Allele1, Allele2, Effect, SE, P, FREQA1) %>%
  rename_with(~ paste0(., "_MA"))

femaleGWAS_GenScot_df_format <- femaleGWAS_GenScot_df %>%
  mutate(MarkerName_match = paste(CHR, POS, A1, A2, sep = ":"),
         MarkerName_match_flipped = paste(CHR, POS, A2, A1, sep = ":")) %>%
  rename(rsID = SNP,
         BP = POS,
         Allele1 = A1,
         Allele2 = A2,
         Effect = BETA,
         FREQA1 = AF1) %>%
  select(MarkerName_match, MarkerName_match_flipped, rsID, CHR, BP, Allele1, Allele2, Effect, SE, P, FREQA1) %>%
  rename_with(~ paste0(., "_GenScot"))



merged_df <- femaleSNPs_MA_df %>%
  inner_join(femaleGWAS_GenScot_df_format, by = c("MarkerName_match_MA" = "MarkerName_match_GenScot")) %>%
  rename(MarkerName = MarkerName_match_MA)

merged_flipped_df <- femaleSNPs_MA_df %>%
  inner_join(femaleGWAS_GenScot_df_format, by = c("MarkerName_match_MA" = "MarkerName_match_flipped_GenScot")) %>%
  rename(MarkerName = MarkerName_match_MA) %>%
  mutate(Effect_GenScot = -1 * Effect_GenScot)

femaleSNPs_MA_GenScot_df <- merged_df %>%
  bind_rows(merged_flipped_df)






### COMPARE SIG SNPs ###

femaleSNPs_MA_GenScot_df <- femaleSNPs_MA_GenScot_df %>%
  mutate(Effect_same_direction = case_when(
    (Effect_MA > 0 & Effect_GenScot > 0) | (Effect_MA < 0 & Effect_GenScot < 0) ~ "Yes",
    (Effect_MA > 0 & Effect_GenScot < 0) | (Effect_MA < 0 & Effect_GenScot > 0) ~ "No"
  ),
  Effect_MA_minus_GenScot = Effect_MA - Effect_GenScot)


femaleSNPs_MA_GenScot_df %>% count(Effect_same_direction)
# 11 in same direction, 5 not
femaleSNPs_MA_GenScot_df %>% count(Effect_MA_minus_GenScot)

outfile <- paste(directory, "Female_ind_lead_sig_SNPs_MA_GenScot.txt", sep="")
write.table(femaleSNPs_MA_GenScot_df, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


### BINOMIAL TEST ###

binom.test(x = 11, n = 16, p = 0.5, alternative = "greater", conf.level = 0.95)











#########MALES##############

### EXTRACT GENOME-WIDE SIG SNPs ###

# Extract only lead independent sig SNPs and match with GenScot sumstats (ensuring alleles matched correctly and if flipped then flip beta)

male_ind_sig_SNPs <- clump_df_male %>%
  filter(P < 5e-08) %>%
  select(SNP) %>%
  pull(SNP)

# Pull out lead ind sig SNPs using MarkerName (same MarkerName used in sumstats and clumping)
# MarkerName is based on dbSNP ref and doesn't necessarily match A1:A2 so make new markername to match with GenScot replication sumstats
maleSNPs_MA_df <- maleGWAS_MA_df %>%
  filter(MarkerName %in% male_ind_sig_SNPs) %>%
  rename(SE = StdErr,
         rsID = rsID_build37,
         FREQA1 = Freq1) %>%
  mutate(MarkerName_match = paste(CHR, BP, toupper(Allele1), toupper(Allele2), sep = ":")) %>%
  select(MarkerName_match, rsID, CHR, BP, Allele1, Allele2, Effect, SE, P, FREQA1) %>%
  rename_with(~ paste0(., "_MA"))

maleGWAS_GenScot_df_format <- maleGWAS_GenScot_df %>%
  mutate(MarkerName_match = paste(CHR, POS, A1, A2, sep = ":"),
         MarkerName_match_flipped = paste(CHR, POS, A2, A1, sep = ":")) %>%
  rename(rsID = SNP,
         BP = POS,
         Allele1 = A1,
         Allele2 = A2,
         Effect = BETA,
         FREQA1 = AF1) %>%
  select(MarkerName_match, MarkerName_match_flipped, rsID, CHR, BP, Allele1, Allele2, Effect, SE, P, FREQA1) %>%
  rename_with(~ paste0(., "_GenScot"))



merged_df <- maleSNPs_MA_df %>%
  inner_join(maleGWAS_GenScot_df_format, by = c("MarkerName_match_MA" = "MarkerName_match_GenScot")) %>%
  rename(MarkerName = MarkerName_match_MA)

merged_flipped_df <- maleSNPs_MA_df %>%
  inner_join(maleGWAS_GenScot_df_format, by = c("MarkerName_match_MA" = "MarkerName_match_flipped_GenScot")) %>%
  rename(MarkerName = MarkerName_match_MA) %>%
  mutate(Effect_GenScot = -1 * Effect_GenScot)

maleSNPs_MA_GenScot_df <- merged_df %>%
  bind_rows(merged_flipped_df)







### COMPARE SIG SNPs ###

maleSNPs_MA_GenScot_df <- maleSNPs_MA_GenScot_df %>%
  mutate(Effect_same_direction = case_when(
    (Effect_MA > 0 & Effect_GenScot > 0) | (Effect_MA < 0 & Effect_GenScot < 0) ~ "Yes",
    (Effect_MA > 0 & Effect_GenScot < 0) | (Effect_MA < 0 & Effect_GenScot > 0) ~ "No"
  ),
  Effect_MA_minus_GenScot = Effect_MA - Effect_GenScot)


maleSNPs_MA_GenScot_df %>% count(Effect_same_direction)
# 5 in same direction, 2 not
maleSNPs_MA_GenScot_df %>% count(Effect_MA_minus_GenScot)

# Sig SNP in chr 23 not in GenScot sumstats - double check this
maleGWAS_GenScot_df_format %>% filter(CHR_GenScot == "23" & BP_GenScot == "28308742")

outfile <- paste(directory, "Male_ind_lead_sig_SNPs_MA_GenScot.txt", sep="")
write.table(maleSNPs_MA_GenScot_df, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)

### BINOMIAL TEST ###

binom.test(x = 5, n = 7, p = 0.5, alternative = "greater", conf.level = 0.95)

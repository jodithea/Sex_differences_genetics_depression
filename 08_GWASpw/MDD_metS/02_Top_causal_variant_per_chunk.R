directory="/path/06_GWAS_pw/MDD_metS/"

femaleMDD_metS_SNP <- "/path/06_GWAS_pw/MDD_metS/femaleMDD_metS/GWAS_pw_femaleMDD_metS.bfs.gz"

maleMDD_metS_SNP <- "/path/06_GWAS_pw/MDD_metS/maleMDD_metS/GWAS_pw_maleMDD_metS.bfs.gz"

femaleMDD_metS_shared_regions <- "/path/06_GWAS_pw/MDD_metS/femaleMDD_metS_model3_vs_maleMDD_metS_model3_regions_female_only.txt"

maleMDD_metS_shared_regions <- "/path/06_GWAS_pw/MDD_metS/femaleMDD_metS_model3_vs_maleMDD_metS_model3_regions_male_only.txt"

bothsexes_MDD_metS_shared_regions <- "/path/06_GWAS_pw/MDD_metS/femaleMDD_metS_model3_vs_maleMDD_metS_model3_regions_both.txt"

### LOAD PACKAGES ###
library(tidyverse)

### READ IN DATA ###

# GWAS-pw data for all SNPs
femaleMDD_metS_SNP_df <- read.table(gzfile(femaleMDD_metS_SNP), header = T, stringsAsFactors = F)
maleMDD_metS_SNP_df <- read.table(gzfile(maleMDD_metS_SNP), header = T, stringsAsFactors = F)

# Regions from gwas-pw identied as support for model 3 (PPA > 0.5) - shared association with both traits
femaleMDD_metS_shared_regions_df <- read.table(femaleMDD_metS_shared_regions, header = T, stringsAsFactors = F)
maleMDD_metS_shared_regions_df <- read.table(maleMDD_metS_shared_regions, header = T, stringsAsFactors = F)
bothsexes_MDD_metS_shared_regions_df <- read.table(bothsexes_MDD_metS_shared_regions, header = T, stringsAsFactors = F)


### TOP CAUSAL VARIANT FOR SHARED REGION MDD - metS: FEMALS ONLY ###

# Within each chunk identified as having a shared causal variant between the 2 traits, identify which is the likely causal variant (SNP with the highest PPA_3 = posterior probability that this SNP is the causal one under model 3)
female_top_SNP_per_chunk <- femaleMDD_metS_SNP_df %>%
  filter(chunk %in% femaleMDD_metS_shared_regions_df$chunk) %>%
  group_by(chunk) %>%
  arrange(-PPA_3) %>%
  slice_head(n = 1) %>%
  ungroup()

outfile <- paste(directory, "Top_SNP_per_chunk_MDD_metS_females_only.txt", sep="")
write.table(female_top_SNP_per_chunk, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)

# Print SNPs per chunk with PPA_3 > 0.5
female_PPA_05_SNP_per_chunk <- femaleMDD_metS_SNP_df %>%
  filter(chunk %in% femaleMDD_metS_shared_regions_df$chunk) %>%
  group_by(chunk) %>%
  filter(PPA_3 > 0.5) %>%
  ungroup()

outfile <- paste(directory, "SNP_with_PPA_above_05_per_chunk_MDD_metS_females_only.txt", sep="")
write.table(female_PPA_05_SNP_per_chunk, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


### TOP CAUSAL VARIANT FOR SHARED REGION MDD - metS: MALES ONLY ###

# Within each chunk identified as having a shared causal variant between the 2 traits, identify which is the likely causal variant (SNP with the highest PPA_3 = posterior probability that this SNP is the causal one under model 3)
male_top_SNP_per_chunk <- maleMDD_metS_SNP_df %>%
  filter(chunk %in% maleMDD_metS_shared_regions_df$chunk) %>%
  group_by(chunk) %>%
  arrange(-PPA_3) %>%
  slice_head(n = 1) %>%
  ungroup()

outfile <- paste(directory, "Top_SNP_per_chunk_MDD_metS_males_only.txt", sep="")
write.table(male_top_SNP_per_chunk, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)

# Print SNPs per chunk with PPA_3 > 0.5
male_PPA_05_SNP_per_chunk <- maleMDD_metS_SNP_df %>%
  filter(chunk %in% maleMDD_metS_shared_regions_df$chunk) %>%
  group_by(chunk) %>%
  filter(PPA_3 > 0.5) %>%
  ungroup()

outfile <- paste(directory, "SNP_with_PPA_above_05_per_chunk_MDD_metS_males_only.txt", sep="")
write.table(male_PPA_05_SNP_per_chunk, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


### TOP CAUSAL VARIANT FOR SHARED REGION MDD - metS: BOTH SEXES ###

bothsexes_top_SNP_per_chunk_female <- femaleMDD_metS_SNP_df %>%
  filter(chunk %in% bothsexes_MDD_metS_shared_regions_df$chunk) %>%
  group_by(chunk) %>%
  arrange(-PPA_3) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename_with(~ paste0(., "_topfemaleSNP"))

bothsexes_top_SNP_per_chunk_male <- maleMDD_metS_SNP_df %>%
  filter(chunk %in% bothsexes_MDD_metS_shared_regions_df$chunk) %>%
  group_by(chunk) %>%
  arrange(-PPA_3) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename_with(~ paste0(., "_topmaleSNP"))


# Female and male top SNP for each chunk can be different so find matching SNP in opposite sex
bothsexes_top_SNP_per_chunk_female <- bothsexes_top_SNP_per_chunk_female %>%
  inner_join(maleMDD_metS_SNP_df, by = join_by("id_topfemaleSNP" == "id"))

outfile <- paste(directory, "Top_SNP_per_chunk_MDD_metS_both_sexes_femaleSNPs.txt", sep="")
write.table(bothsexes_top_SNP_per_chunk_female, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



bothsexes_top_SNP_per_chunk_male <- bothsexes_top_SNP_per_chunk_male %>%
  inner_join(femaleMDD_metS_SNP_df, by = join_by("id_topmaleSNP" == "id"))

outfile <- paste(directory, "Top_SNP_per_chunk_MDD_metS_both_sexes_maleSNPs.txt", sep="")
write.table(bothsexes_top_SNP_per_chunk_male, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)



# Print SNPs per chunk with PPA_3 > 0.5 in both females and males

bothsexes_PPA_05_SNP_per_chunk_female <- femaleMDD_metS_SNP_df %>%
  filter(chunk %in% bothsexes_MDD_metS_shared_regions_df$chunk) %>%
  group_by(chunk) %>%
  filter(PPA_3 > 0.5) %>%
  ungroup() %>%
  rename_with(~ paste0(., "_female"))

bothsexes_PPA_05_SNP_per_chunk_male <- maleMDD_metS_SNP_df %>%
  filter(chunk %in% bothsexes_MDD_metS_shared_regions_df$chunk) %>%
  group_by(chunk) %>%
  filter(PPA_3 > 0.5) %>%
  ungroup() %>%
  rename_with(~ paste0(., "_male"))

bothsexes_PPA_05_SNP_per_chunk <- bothsexes_PPA_05_SNP_per_chunk_female %>%
  inner_join(bothsexes_PPA_05_SNP_per_chunk_male, by = join_by("id_female" == "id_male")) %>%
  rename("id" = "id_female")

outfile <- paste(directory, "SNP_with_PPA_above_05_per_chunk_MDD_metS_both_sexes.txt", sep="")
write.table(bothsexes_PPA_05_SNP_per_chunk, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)


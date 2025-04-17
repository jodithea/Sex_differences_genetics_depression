# Plot mean of SE vs effective sample size of each cohort - expect SE to reflect Neff

directory = "/path/01_Format_sumstats/Males/"

GWAS_directory = "/path/01_Format_sumstats/"

AGDS_male = paste0(GWAS_directory, "Males/AGDS_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

AllOfUs_male = paste0(GWAS_directory, "Males/AllOfUs_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

Bionic_male = paste0(GWAS_directory, "Males/Bionic_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

Blokland_male = paste0(GWAS_directory, "Males/Blokland21_MDD_male_sumstats_formatted_formetaanalysis.txt")

GLAD_male = paste0(GWAS_directory, "Males/GLAD_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

UKB_male = paste0(GWAS_directory, "Males/UKB_MDD_male_sumstats_formatted_formetaanalysis.txt", sep = "")

## ---- LoadPackages
library(tidyverse)
## ----

## ---- LoadData
load_dataframe <- function(df_name) {
  df <- read.table(df_name, header = T, stringsAsFactors = F)
  
  return(df)
}

file_paths <- c(
  AGDS_male = AGDS_male,
  AllOfUs_male = AllOfUs_male,
  Bionic_male = Bionic_male,
  Blokland_male = Blokland_male,
  GLAD_male = GLAD_male,
  UKB_male = UKB_male
)

data_frames <- lapply(file_paths, load_dataframe)
# Access individual data frames
# data_frames[["AGDS_male"]]
## ----

## ---- Dataframe
# Calculate Neff and mean of beta SE for each cohort using just SNPs with freq of 0.45 - 0.55

process_dataset <- function(dataset, cohort_name, N, P) {
  dataset %>%
    summarise(Cohort = cohort_name, 
              Mean_SE_0_0.2 = mean(SE[FREQA1 > 0 & FREQA1 < 0.2], na.rm = TRUE),
              Mean_SE_0.2_0.45 = mean(SE[FREQA1 > 0.2 & FREQA1 < 0.45], na.rm = TRUE),
              Mean_SE_0.45_0.55 = mean(SE[FREQA1 > 0.45 & FREQA1 < 0.55], na.rm = TRUE),
              Mean_SE_all = mean(SE, na.rm = TRUE)) %>%
    mutate(Neff = 4 * N * P * (1 - P)) %>%
    pivot_longer(cols = starts_with("Mean_SE"), 
                 names_to = "FREQA1_Bracket", 
                 values_to = "Mean_SE") %>%
    mutate(FREQA1_Bracket = case_when(
      FREQA1_Bracket == "Mean_SE_0_0.2" ~ "0-0.2",
      FREQA1_Bracket == "Mean_SE_0.2_0.45" ~ "0.2-0.45",
      FREQA1_Bracket == "Mean_SE_0.45_0.55" ~ "0.45-0.55",
      FREQA1_Bracket == "Mean_SE_all" ~ "All"
    ))
}

# Define cohort name, N and P for each dataset/cohort
cohort_details <- list(
  AGDS = list(dataset = data_frames[["AGDS_male"]], N = 9775, P = 0.3247),
  AllOfUs = list(dataset = data_frames[["AllOfUs_male"]], N = 46972, P = 0.2371),
  Bionic = list(dataset = data_frames[["Bionic_male"]], N = 24445, P = 0.1813),
  Blokland = list(dataset = data_frames[["Blokland_male"]], N = 28993, P = 0.6484),
  GLAD = list(dataset = data_frames[["GLAD_male"]], N = 7681, P = 0.6062),
  UKB = list(dataset = data_frames[["UKB_male"]], N = 79124, P = 0.2857)
)


result <- lapply(names(cohort_details), function(cohort_name) {
  cohort <- cohort_details[[cohort_name]] # Extract details for the current cohort
  process_dataset(
    dataset = cohort$dataset,
    cohort_name = cohort_name,
    N = cohort$N,
    P = cohort$P
  )
})


# Combine the results into a single dataframe
summary_df <- do.call(rbind, result)
## ----

## ---- Plot
plot <- ggplot(summary_df, aes(x = Neff, y = Mean_SE)) +
  geom_point() +
  geom_text(aes(label = Cohort), nudge_y = 0.002) + 
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ FREQA1_Bracket, scales = "free") +
  scale_x_continuous(limits = c(0, 70000),
                     breaks = seq(0, 70000, 10000),
                     labels = seq(0, 70000, 10000)) +
  ggtitle("Mean SE vs Neff for SNPs by FREQA1 Bracket") +
  theme_classic()

plot

outfile = paste0(directory, "/SE_vs_Neff_male_cohorts.png")
ggsave(plot, file = outfile, width = 30, height = 30, unit = "cm")
## ----

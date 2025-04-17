# This script compare h2 of females vs males in the UK Biobank full sumstats
# Using lifetime pop prev of 20% for females and 10% for males and h2 calculated in SBayesS

directory = "/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_full/"

UKB_full_MCMC = "/path/12_SBayesS_h2/UKB_all_MCMC_female_male.RData"


## PACKAGES ##
library(coda)       # HPDI
library(tidyverse)



## LOAD DATA ##

UKB_full_MCMC <- load(file = UKB_full_MCMC)




### Liability h2 for a range of K ###
# Create a list of column names, sex labels, and K values for each column
liability_columns <- c("hsq_F_L25", "hsq_F_L20", "hsq_F_L15", "hsq_F_L10", "hsq_F_L5",
                       "hsq_M_L25", "hsq_M_L20", "hsq_M_L15", "hsq_M_L10", "hsq_M_L5")

sex_labels <- c(rep("female", 5), rep("male", 5))
K_values <- c(0.25, 0.20, 0.15, 0.10, 0.05, 0.25, 0.20, 0.15, 0.10, 0.05)


# Initialize an empty list to store results for liability scale
results_list <- list()

# Iterate over each column and calculate posterior mean and 95% HPD interval for liability scale
for (i in seq_along(liability_columns)) {
  column_name <- liability_columns[i]
  sex <- sex_labels[i]
  K <- K_values[i]
  
  # Calculate the posterior mean for liability scale
  # 2,000 data points = chain length of 25,000 with thinning of 10 and burnin of 5,000 already removed
  posterior_mean_liability <- UKB_all_mcmc %>%
    summarize(mean = mean(.data[[column_name]])) %>%
    pull(mean)
  
  # Convert to MCMC object and calculate 95% HPD interval for liability scale
  mcmc_obj_liability <- as.mcmc(UKB_all_mcmc %>% pull(.data[[column_name]]))
  hpd_interval_liability <- HPDinterval(mcmc_obj_liability, prob = 0.95)[1, ]
  
  # Store liability scale results in the list
  results_list[[i]] <- data.frame(
    mean_hsq = posterior_mean_liability,
    hpdi_lower = hpd_interval_liability[1],
    hpdi_upper = hpd_interval_liability[2],
    sex = sex,
    K = K
  )
}

# Combine all liability-scale results into a single data frame
h2_liability_varying_K <- do.call(rbind, results_list)

### Add the observed-scale heritability for female and male ###
# Female observed scale
observed_female <- UKB_all_mcmc %>%
  summarize(mean = mean(hsq_female)) %>%
  pull(mean)

observed_female_hpd <- HPDinterval(as.mcmc(UKB_all_mcmc %>% pull(hsq_female)), prob = 0.95)[1, ]

# Male observed scale
observed_male <- UKB_all_mcmc %>%
  summarize(mean = mean(hsq_male)) %>%
  pull(mean)

observed_male_hpd <- HPDinterval(as.mcmc(UKB_all_mcmc %>% pull(hsq_male)), prob = 0.95)[1, ]

# Create data frame for observed values and bind to liability scale results
observed_h2 <- data.frame(
  mean_hsq = c(observed_female, observed_male),
  hpdi_lower = c(observed_female_hpd[1], observed_male_hpd[1]),
  hpdi_upper = c(observed_female_hpd[2], observed_male_hpd[2]),
  sex = c("female", "male"),
  K = "observed"
)


# Save the combined data frame
save(h2_liability_varying_K, file = paste0(directory, "UKB_full_sumstats_h2_liability_varying_K.RData"))
save(observed_h2, file = paste0(directory, "UKB_full_sumstats_observed_h2.RData"))

write.table(h2_liability_varying_K, file = paste0(directory, "UKB_full_sumstats_h2_liability_varying_K.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(observed_h2, file = paste0(directory, "UKB_full_sumstats_observed_h2.txt"), col.names = T, row.names = F, sep = "\t", quote = F)




## POSTERIOR PROB FEMALE > MALE ## 

# Count the frequency in the MCMC sample that h2_female > h2_male = posterior probability that female h2 is larger than male h2
post_prob_M_F <- UKB_all_mcmc %>%
  mutate(hsq_f_minus_m = hsq_female - hsq_male) %>%
  summarize(Perc = (sum(hsq_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

post_prob_M_F_liability <- UKB_all_mcmc %>%
  mutate(hsq_f_minus_m = hsq_F_L20 - hsq_M_L10) %>%
  summarize(Perc = (sum(hsq_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)

# Print the result
sink(paste0(directory, "/male_vs_female_h2_UKB_full_sumstats.txt"))
cat("Posterior probability that h2_female (observed scale) > h2_male (observed scale):", post_prob_M_F, "\n")
cat("Posterior probability that h2_female (liability scale with K = 0.20) > h2_male (liability scale with K = 0.10):", post_prob_M_F_liability, "\n")
sink()

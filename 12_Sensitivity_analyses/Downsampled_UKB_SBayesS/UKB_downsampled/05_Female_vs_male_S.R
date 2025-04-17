# This script compare selection parameter of females vs males in the UK Biobank full sumstats


directory = "/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/"

UKB_downsampled_MCMC <- paste0(directory, "UKB_downsampled_all_MCMC_female_male.RData", sep = "")



## PACKAGES ##
library(coda)       # HPDI
library(tidyverse)


## LOAD DATA ##

UKB_downsampled_MCMC <- load(file = UKB_downsampled_MCMC)



### Selection parameter for female and male ###
# Female
female <- UKB_all_mcmc %>%
  summarize(mean = mean(S_female)) %>%
  pull(mean)

female_hpd <- HPDinterval(as.mcmc(UKB_all_mcmc %>% pull(S_female)), prob = 0.95)[1, ]

# Male
male <- UKB_all_mcmc %>%
  summarize(mean = mean(S_male)) %>%
  pull(mean)

male_hpd <- HPDinterval(as.mcmc(UKB_all_mcmc %>% pull(S_male)), prob = 0.95)[1, ]

# Create data frame for values
S <- data.frame(
  mean_S = c(female, male),
  hpdi_lower = c(female_hpd[1], male_hpd[1]),
  hpdi_upper = c(female_hpd[2], male_hpd[2]),
  sex = c("female", "male")
)


# Save
save(S, file = paste0(directory, "UKB_downsampled_sumstats_selection.RData"))

write.table(S, file = paste0(directory, "UKB_downsampled_sumstats_selection.txt"), col.names = T, row.names = F, sep = "\t", quote = F)





## POSTERIOR PROB FEMALE > MALE ##

# Count the frequency in the MCMC sample that S_female > S_male = posterior probability that female S is larger than male S
post_prob_M_F <- UKB_all_mcmc %>%
  mutate(S_f_minus_m = S_female - S_male) %>%
  summarize(Perc = (sum(S_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)


# Print the result
sink(paste0(directory, "/male_vs_female_S_UKB_downsampled_sumstats.txt"))
cat("Posterior probability that S_female > S_male:", post_prob_M_F, "\n")
sink()

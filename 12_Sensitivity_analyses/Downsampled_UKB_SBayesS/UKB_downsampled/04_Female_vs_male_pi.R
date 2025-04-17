# This script compare polygenicity of females vs males in the UK Biobank downsampled sumstats

directory = directory = "/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/"

UKB_downsampled_MCMC <- paste0(directory, "UKB_downsampled_all_MCMC_female_male.RData", sep = "")


## PACKAGES ##
library(coda)       # HPDI
library(tidyverse)


## LOAD DATA ##

UKB_downsampled_MCMC <- load(file = UKB_downsampled_MCMC)



### Polygenicity for female and male ###
# Female
female <- UKB_all_mcmc %>%
  summarize(mean = mean(Pi_female)) %>%
  pull(mean)

female_hpd <- HPDinterval(as.mcmc(UKB_all_mcmc %>% pull(Pi_female)), prob = 0.95)[1, ]

# Male
male <- UKB_all_mcmc %>%
  summarize(mean = mean(Pi_male)) %>%
  pull(mean)

male_hpd <- HPDinterval(as.mcmc(UKB_all_mcmc %>% pull(Pi_male)), prob = 0.95)[1, ]

# Create data frame for values
pi <- data.frame(
  mean_pi = c(female, male),
  hpdi_lower = c(female_hpd[1], male_hpd[1]),
  hpdi_upper = c(female_hpd[2], male_hpd[2]),
  sex = c("female", "male")
)


# Save
save(pi, file = paste0(directory, "UKB_downsampled_sumstats_polygenicity.RData"))

write.table(pi, file = paste0(directory, "UKB_downsampled_sumstats_polygenicity.txt"), col.names = T, row.names = F, sep = "\t", quote = F)





## POSTERIOR PROB FEMALE > MALE ##

# Count the frequency in the MCMC sample that pi_female > pi_male = posterior probability that female pi is larger than male pi
post_prob_M_F <- UKB_all_mcmc %>%
  mutate(pi_f_minus_m = Pi_female - Pi_male) %>%
  summarize(Perc = (sum(pi_f_minus_m > 0)/n())*100) %>%
  mutate(Perc = round(Perc, 2)) %>%
  pull(Perc)


# Print the result
sink(paste0(directory, "/male_vs_female_pi_UKB_downsampled_sumstats.txt"))
cat("Posterior probability that pi_female > pi_male:", post_prob_M_F, "\n")
sink()

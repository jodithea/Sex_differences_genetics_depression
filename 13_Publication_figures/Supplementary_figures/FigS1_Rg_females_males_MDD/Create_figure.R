directory = "/path/10_Publication_figures/Supplementary_figures/Rg_females_males_MDD/"

rg_female="/path/04_LDSC/SNPrg/Females/rg_femaleMDD_correlation_table.txt"

rg_male="/path/04_LDSC/SNPrg/Males/rg_maleMDD_correlation_table.txt"

### LOAD PACKAGES ###

library(patchwork)
library(tidyverse)



### LOAD DATA ###

rg_female_df <- read.table(rg_female, header = T, stringsAsFactors = F)

rg_male_df <- read.table(rg_male, header = T, stringsAsFactors = F)

rg_df <- bind_rows(rg_female_df, rg_male_df) %>%
  mutate(CI = 1.96 * se,
         p1 = gsub("_sumstats", "", p1),
         p2 = gsub("_sumstats", "", p2))


### RG WITH PREVIOUS MDD ###

order <- c('AdamsMDD',
           'BloklandMDDfemale',
           'BloklandMDDmale',
           'SilveiraMDDfemale',
           'SilveiraMDDmale')
order <- rev(order)

plot_df <- rg_df %>%
  filter(p2 == 'AdamsMDD' |
           p2 == 'BloklandMDDfemale' |
           p2 == 'BloklandMDDmale' |
           p2 == 'SilveiraMDDfemale' |
           p2 == 'SilveiraMDDmale') %>%
  mutate(p2 = factor(p2, levels = order),
         Sex = factor(p1, levels = c("Allcohorts_female_MDD", "Allcohorts_male_MDD")))



rg_plot_MDD <- ggplot(plot_df, aes(x = rg, y = p2, colour = p1)) +
  geom_pointrange(aes(x = rg, xmin = rg + CI, xmax = rg - CI),
                  position = position_dodge(width = -0.8)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  scale_colour_viridis_d("Sex", begin = 1, end = 0,
                         labels = c("Females", "Males")) +
  scale_y_discrete(breaks = c('AdamsMDD',
                              'BloklandMDDfemale',
                              'BloklandMDDmale',
                              'SilveiraMDDfemale',
                              'SilveiraMDDmale'),
                   labels = c('Adams et al., 2025',
                              'Blokland et al., 2022: Females',
                              'Blokland et al., 2022: Males',
                              'Silveira et al., 2023: Females',
                              'Silveira et al., 2023: Males')) +
  labs(x = "Genetic Correlation",
       y = "") +
  annotate("text", x = 1.1, y = "AdamsMDD", label = "*", size = 5) +
  annotate("text", x = 1.3, y = "BloklandMDDfemale", label = "*", size = 5) +
  annotate("text", x = 1.7, y = "BloklandMDDmale", label = "*", size = 5) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        legend.position = "top",
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))


rg_plot_MDD

outfile <- paste(directory, "Rg_female_vs_male_MDD.png", sep="")
ggsave(rg_plot_MDD, width = 20, height = 10, unit = "cm", file = outfile)

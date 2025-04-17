directory = "/path/10_Publication_figures/Supplementary_figures/GxS_forest_plot/"

GxS="/path/03_Metal/GxS/fulldc/Metaanalysis_MDD_GxS_fulldc_AllCohorts_QCed_rsID.txt"

clumpfile_inter = "/path/05_Clumping/GxS/fulldc/clumped_GxS_fulldc_all_results.txt"

femaleGWAS = "/path/03_Metal/Females/Metaanalysis_MDD_female_AllCohorts_QCed_rsID.txt"

maleGWAS = "/path/03_Metal/Males/Metaanalysis_MDD_male_AllCohorts_QCed_rsID.txt"

clumpfile_female = "/path/05_Clumping/Females/clumped_female_all_results.txt"

clumpfile_male = "/path/05_Clumping/Males/clumped_male_all_results.txt"


### Packages ###
library(patchwork)
library(tidyverse)


### Load Data ###
GxS_df <- read.table(GxS, header = TRUE, stringsAsFactors = F)

clump_df_inter <- read.table(clumpfile_inter, header = TRUE, stringsAsFactors = F

femaleGWAS_df <- read.table(femaleGWAS, header = T, stringsAsFactors = F)

maleGWAS_df <- read.table(maleGWAS, header = T, stringsAsFactors = F)

clump_df_female <- read.table(clumpfile_female, header = TRUE, stringsAsFactors = F)

clump_df_male <- read.table(clumpfile_male, header = TRUE, stringsAsFactors = F)


### FOREST PLOT ###

# Independent lead nominally sig SNPs from GxS analysis

clumpinter_sig <- clump_df_inter %>%
  filter(P < 1e-06) %>%
  select(SNP) %>%
  rename(MarkerName = SNP)


# Extract these SNPs for both female and male sumstas
femaleGWAS_clumpinter_sig <- femaleGWAS_df %>%
  inner_join(clumpinter_sig) %>%
  mutate(Dataset = "female")

maleGWAS_df_clumpinter_sig <- maleGWAS_df %>%
  inner_join(clumpinter_sig) %>%
  mutate(Dataset = "male")


# Create one df with all data to make plot, ordered by effect size in female sumstas
df_plot <- femaleGWAS_clumpinter_sig %>%
  full_join(maleGWAS_df_clumpinter_sig)

df_plot <- df_plot %>%
  arrange(CHR) %>%
  mutate(MarkerName = factor(MarkerName, levels = unique(MarkerName)),
         CI = 1.96 * StdErr)


forest_plot <- ggplot(df_plot, aes(x = Effect, y = rsID_build37, colour = Dataset)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dotted') +
  geom_linerange(aes(xmin = Effect - CI, xmax = Effect + CI),
                 position = position_dodge(width = -0.8)) +
  geom_point(position = position_dodge(width = -0.8)) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Females", "Males")) +
  scale_x_continuous("Effect Size",
                     lim = c(-0.10, 0.1),
                     breaks = seq(-0.1, 0.1, 0.05),
                     labels = seq(-0.1, 0.1, 0.05)) +
  scale_y_discrete("",
                   breaks = c('rs6080675',
                              'rs28573687',
                              'rs12312238',
                              'rs12092435'),
                   limits = c('rs6080675',
                              'rs28573687',
                              'rs12312238',
                              'rs12092435')) +
  theme_classic() +
  theme(legend.position = "top",
        legend.key.spacing.x = unit(0.5, 'cm'),
        legend.title = element_blank(),
        text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(0,10,0,0), "pt"))

forest_plot

outfile <- paste(directory, "GxS_nominally_sig_SNPs_forest_plot.png", sep="")
ggsave(forest_plot, width = 21, height = 10, unit = "cm", file = outfile)


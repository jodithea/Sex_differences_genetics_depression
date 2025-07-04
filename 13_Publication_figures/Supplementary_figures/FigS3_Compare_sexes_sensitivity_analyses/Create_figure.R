# Create publication figure of sensitivity analyses comparing males and females
# SBayesS (heritability, polygenicity, selection paramter) downsampled and across-cohort heterogeneity, h2 varying population prevalence, h2 with varying unscreened controls in males, mixer with downsampling

directory = "/path/10_Publication_figures/Supplementary_figures/Compare_sexes_sensitivity_analyses/"



### Packages ###
library(patchwork)
library(ggforce)      # ggplot circles
library(tidyverse)

# heritability data
load("/path/12_SBayesS_h2/GWAS_MA_sumstats_h2_liability_varying_K.RData")
h2_liability_varying_K_MA <- h2_liability_varying_K

load("/path/12_SBayesS_h2/GWAS_MA_sumstats_h2_liability_varying_u_males.RData")

load("/path/12_SBayesS_h2/GWAS_MA_sumstats_observed_h2.RData")
observed_h2_MA <- observed_h2

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_full/UKB_full_sumstats_h2_liability_varying_K.RData")
h2_liability_varying_K_UKB_full <- h2_liability_varying_K

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/UKB_downsampled_sumstats_h2_liability_varying_K.RData")
h2_liability_varying_K_UKB_downsampled <- h2_liability_varying_K


load("/path/12_SBayesS_h2/Cohorts_h2_liability_varying_K.RData")

# polygneicity data

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_full/UKB_full_sumstats_polygenicity.RData")
pi_UKB_full <- pi

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/UKB_downsampled_sumstats_polygenicity.RData")
pi_UKB_downsampled <- pi

load("/path/12_SBayesS_h2/Pi_data_MA_cohorts.RData")
pi_data_cohorts <- all_data_pi %>% 
  filter(Type == "cohorts")

# Selection parameter data

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_full/UKB_full_sumstats_selection.RData")
S_UKB_full <- S

load("/path/08_Sensitivity_analyses/Downsample/SBayesS/UKB_downsampled/UKB_downsampled_sumstats_selection.RData")
S_UKB_downsampled <- S

load("/path/12_SBayesS_h2/S_data_MA_cohorts.RData")
S_data_cohorts <- all_data_S %>% 
  filter(Type == "cohorts")



### SNP-based heritability Plot: Liability scale with UKB full and downsampled ###

h2_liability_UKB_full <- h2_liability_varying_K_UKB_full %>%
  filter((sex == "female" & K == 0.20) |
           (sex == "male" & K == 0.10)) %>%
  mutate(cohort = "UKB_full") %>%
  select(-K)

h2_liability_UKB_downsampled <- h2_liability_varying_K_UKB_downsampled %>%
  filter((sex == "female" & K == 0.20) |
           (sex == "male" & K == 0.10)) %>%
  mutate(cohort = "UKB_downsampled") %>%
  select(-K)

h2_liability_UKB <- h2_liability_UKB_full %>%
  bind_rows(h2_liability_UKB_downsampled)

h2_plot_sexes_UKB <- ggplot(h2_liability_UKB, aes(x = cohort, y = mean_hsq * 100, colour = sex)) +
  geom_point(position = position_dodge(width = -0.4)) +
  geom_linerange(aes(ymin = hpdi_lower * 100, ymax = hpdi_upper * 100),
                 position = position_dodge(width = -0.4)) +
  scale_y_continuous(bquote('Liability ' ~{h^2} [SNP]~' (%)'),
                     limits = c(10, 22),
                     breaks = seq(10, 20, 2)) +
  scale_x_discrete("",
                   limits = c("UKB_full", "UKB_downsampled"),
                   labels = c("UK B \nfull", "UK B \ndownsampled")) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Female", "Male")) +
  guides(color = "none") +
  annotate("text", x = 1, y = 21, label = "P(F > M) \n= 99%", size = 3) +
  annotate("text", x = 2, y = 21, label = "P(F > M) \n= 100%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))
h2_plot_sexes_UKB


### Polygenicity Plot: UKB full and downsampled ###

pi_UKB <- (pi_UKB_full %>% mutate(cohort = "UKB_full")) %>%
  bind_rows(pi_UKB_downsampled %>% mutate(cohort = "UKB_downsampled"))


pi_plot_sexes_UKB <- ggplot(pi_UKB, aes(x = cohort, y = mean_pi, colour = sex)) +
  geom_point(position = position_dodge(width = -0.4)) +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper),
                 position = position_dodge(width = -0.4)) +
  scale_y_continuous("Polygenicity (\u03c0)",
                     limits = c(0, 0.024),
                     breaks = seq(0, 0.02, 0.004)) +
  scale_x_discrete("",
                   limits = c("UKB_full", "UKB_downsampled"),
                   labels = c("UK B \nfull", "UK B \ndownsampled")) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Female", "Male")) +
  guides(color = "none") +
  annotate("text", x = 1, y = 0.022, label = "P(F > M) \n= 87%", size = 3) +
  annotate("text", x = 2, y = 0.022, label = "P(F > M) \n= 29%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))
pi_plot_sexes_UKB


### Selection Parameter Plot: UKB full and downsampled ###

S_UKB <- (S_UKB_full %>% mutate(cohort = "UKB_full")) %>%
  bind_rows(S_UKB_downsampled %>% mutate(cohort = "UKB_downsampled"))


S_plot_sexes_UKB <- ggplot(S_UKB, aes(x = cohort, y = mean_S, colour = sex)) +
  geom_point(position = position_dodge(width = -0.4)) +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper),
                 position = position_dodge(width = -0.4)) +
  scale_y_continuous(expression("Selection parameter (" * italic(S) * ")"),
                     limits = c(-0.85, 0.4),
                     breaks = seq(-0.8, 0.2, 0.2)) +
  scale_x_discrete("",
                   limits = c("UKB_full", "UKB_downsampled"),
                   labels = c("UK B \nfull", "UK B \ndownsampled")) +
  scale_colour_manual(values = c('#FDE725FF', '#440154FF'),
                      labels = c("Female", "Male")) +
  guides(color = "none") +
  annotate("text", x = 1, y = 0.3, label = "P(F > M) \n= 79%", size = 3) +
  annotate("text", x = 2, y = 0.3, label = "P(F > M) \n= 25%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))
S_plot_sexes_UKB




### SNP-based heritability Plot: Across cohort heterogeneity ###

h2_liability <- h2_liability_varying_K_cohorts %>%
  filter((sex == "female" & K == 0.20) |
           (sex == "male" & K == 0.10)) %>%
  mutate(scale = "Liability") %>%
  select(-K)


h2_plot_sexes_cohorts <- ggplot(h2_liability, aes(x = sex, y = mean_hsq * 100, colour = sex)) +
  geom_point() +
  geom_linerange(aes(ymin = hpdi_lower * 100, ymax = hpdi_upper * 100)) +
  scale_y_continuous(bquote('Liability ' ~{h^2} [SNP]~' (%)'),
                     limits = c(0, 50)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_manual(
    values = c('female' = '#FDE725FF', 
               'male' = '#440154FF'),
    breaks = c("male", "female"),
    labels = c('male' = "Male", 
               'female' = "Female")
  ) +
  annotate("text", x = 1.5, y = 50, label = "P(F > M) = 93%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_plot_sexes_cohorts


### Polygenicity Plot: Across-cohort heterogeneity ###


pi_plot_sexes_cohorts <- ggplot(pi_data_cohorts, aes(x = sex, y = mean_pi, colour = sex)) +
  geom_point(position = position_dodge(width = -0.8)) +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper),
                 position = position_dodge(width = -0.8)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_viridis_d(end = 0, begin = 1) +
  guides(color = "none") +
  scale_y_continuous("Polygenicity (\u03c0)",
                     limits = c(0,0.031)) +
  annotate("text", x = 1.5, y = 0.03, label = "P(F > M) = 79%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

pi_plot_sexes_cohorts

### Selection Plot: Across-cohort heterogeneity ###

S_plot_sexes_cohorts <- ggplot(S_data_cohorts, aes(x = sex, y = mean_S, colour = sex)) +
  geom_point(position = position_dodge(width = -0.8)) +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper),
                 position = position_dodge(width = -0.8)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_viridis_d(end = 0, begin = 1) +
  guides(color = "none") +
  scale_y_continuous(expression("Selection parameter (" * italic(S) * ")"),
                     limits = c(-1.5, 2.5),
                     breaks = seq(-1.5, 2.5, 0.5)) +
  annotate("text", x = 1.5, y = 2.5, label = "P(F > M) = 24%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

S_plot_sexes_cohorts


### SNP-based heritability Plot: Multiple K (SBayesS) ###

# Using results from SBayesS on GWAS sex stratified meta-analysis results

h2_liability_varying_K_MA_plot <- h2_liability_varying_K_MA %>%
  mutate(sex = factor(sex, levels = c("male", "female")))

h2_varying_K_plot <- ggplot(h2_liability_varying_K_MA_plot, aes(x = K, y = mean_hsq*100, colour = sex)) +
  geom_vline(xintercept = 0.10,
             colour = "grey30",
             linetype = 'dashed') +
  geom_vline(xintercept = 0.20,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point(position = position_dodge(width = 0.01)) +
  geom_linerange(aes(ymin = hpdi_lower*100, ymax = hpdi_upper*100),
                 position = position_dodge(width = 0.01)) +
  scale_colour_viridis_d(begin = 0, end = 1,
                         labels = c("Male", "Female")) +
  guides(color = "none") +
  scale_x_continuous("Population Prevalence",
                     breaks = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     labels = c(0.05, 0.10, 0.15, 0.20, 0.25)) +
  scale_y_continuous(bquote('Liability ' ~{h^2} [SNP]~' (%)'),
                     limits = c(6, 14),
                     breaks = seq(6, 14, 2)) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_varying_K_plot


### SNP-based heritability Plot: Varying u (SBayesS) ###

h2_liability_varying_u <- h2_liability_varying_u_males %>%
  mutate(sex = "male") %>%
  slice_tail(n = 11) %>%
  bind_rows(h2_liability_varying_K_MA %>% filter(sex == "female" & K == 0.20) %>% mutate(u = 0))

h2_varying_u_plot <- ggplot(h2_liability_varying_u, aes(x = u, y = mean_hsq*100, colour = sex)) +
  geom_point(position = position_dodge(width = 0.01)) +
  geom_linerange(aes(ymin = hpdi_lower*100, ymax = hpdi_upper*100),
                 position = position_dodge(width = 0.01)) +
  scale_colour_viridis_d(begin = 1, end = 0) +
  guides(color = "none") +
  scale_x_continuous("Proportion of unscreened control subjects",
                     breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(expression("Liability " * italic(h)^2 * " (%)"),
                     limits = c(8, 20),
                     breaks = seq(8, 20, 2)) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))

h2_varying_u_plot



### Mixer Venn: Full UKB ###

# Data: sizes of each group
one <- 3725  # female exclusive
two <- 104   # male exclusive
overlap <- 7516  # Shared between female and male

# Calculate circle areas
area_one <- one + overlap
area_two <- two + overlap

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
rm(pi)
r_1 <- sqrt(area_one / pi)
r_2 <- sqrt(area_two / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(area_one, area_two)
d <- (r_1 + r_2) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_1, r_2),
  group = c("MDD", "BMI")
)

# Plot the Euler diagram with numbers positioned correctly
mixer_venn_UKB_full <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
  guides(fill = "none") +
  coord_fixed(expand = TRUE) +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_1 / 2,
           y = circle_data$y[1],
           label = one,
           size = 3, color = "black", family = "Calibri") +
  annotate("text",
           x = circle_data$x[2] + r_2 / 2,
           y = circle_data$y[2],
           label = two,
           size = 3, color = "black", family = "Calibri") +
  annotate("text",
           x = mean(circle_data$x) + 0.6,
           y = circle_data$y[1],
           label = overlap,
           size = 3, color = "black", family = "Calibri") +
  ggtitle("UK B Full") +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")


mixer_venn_UKB_full



### Mixer Venn: Downsampled UKB ###

# Data: sizes of each group
one <- 5769  # female exclusive
two <- 132   # male exclusive
overlap <- 6645  # Shared between female and male

# Calculate circle areas
area_one <- one + overlap
area_two <- two + overlap

# Convert areas to radii (since area = π * r^2, radius = sqrt(area / π))
rm(pi)
r_1 <- sqrt(area_one / pi)
r_2 <- sqrt(area_two / pi)

# Calculate the proportional overlap distance
overlap_ratio <- overlap / min(area_one, area_two)
d <- (r_1 + r_2) * (1 - overlap_ratio / 2)  # Adjust overlap distance

# Create circle data with calculated positions
circle_data <- data.frame(
  x = c(0, d),  # x-coordinates for each circle
  y = c(0, 0),  # Keep y-coordinates the same for both circles
  radius = c(r_1, r_2),
  group = c("MDD", "BMI")
)

# Plot the Euler diagram with numbers positioned correctly
mixer_venn_UKB_down <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF")) +
  guides(fill = "none") +
  coord_fixed(expand = TRUE) +
  # Add labels inside each circle and in the overlapping section
  annotate("text",
           x = circle_data$x[1] - r_1 / 2,
           y = circle_data$y[1],
           label = one,
           size = 3, color = "black", family = "Calibri") +
  annotate("text",
           x = circle_data$x[2] + r_2 / 2,
           y = circle_data$y[2],
           label = two,
           size = 3, color = "black", family = "Calibri") +
  annotate("text",
           x = mean(circle_data$x) + 0.6,
           y = circle_data$y[1],
           label = overlap,
           size = 3, color = "black", family = "Calibri") +
  ggtitle("UK B Downsampled") +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")


mixer_venn_UKB_down


### Create Panelled Figure ###

figure <- (guide_area() /
             ((h2_plot_sexes_UKB | pi_plot_sexes_UKB | S_plot_sexes_UKB ) + plot_layout(widths = c(1, 1, 1))) /
             ((h2_plot_sexes_cohorts | pi_plot_sexes_cohorts | S_plot_sexes_cohorts) + plot_layout(widths = c(1, 1, 1))) /
             ((h2_varying_K_plot | h2_varying_u_plot) + plot_layout(widths = c(1, 1.5))) /
             ((mixer_venn_UKB_full | mixer_venn_UKB_down))) +
  plot_layout(heights = c(0.2, 1, 1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E", "F", "G", "H", "I", ""))) & 
  theme(legend.position = "top")

figure


outfile <- paste(directory, "Compare_sexes_sensitivity_analyses.png", sep="")
ggsave(figure, width = 21, height = 29, unit = "cm", file = outfile)





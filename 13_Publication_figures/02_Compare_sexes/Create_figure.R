# Create publication figure of tests comparing males and females
# SBayesS (heritability, polygenicity, selection paramter), rg, correlation, polygenic overlap (rg, mixer)

directory = "/path/10_Publication_figures/02_Compare_sexes/"



### Packages ###
library(patchwork)
library(ggforce)      # ggplot circles
library(tidyverse)


### Load Data ###

# heritability data
load("/path/12_SBayesS_h2/GWAS_MA_sumstats_h2_liability_varying_K.RData")
h2_liability_varying_K_MA <- h2_liability_varying_K


# Polygenicity data
load("/path/12_SBayesS_h2/Pi_data_MA_cohorts.RData")


# Selection parameter data
load("/path/12_SBayesS_h2/S_data_MA_cohorts.RData")


# Genetic correlation data
load("/path/04_LDSC/SNPrg/Female_vs_male_all_cohorts/all_data_for_plot.RData")
all_data_rg <- all_data

# Pearson correlation data
load("/path/09_Effect_sizes_plots/Adams_SNPs/all_data_for_plot.RData")
all_data_R <- all_data


### SNP-based heritability Plot: Liability scale ###

# Using results from SBayesS on GWAS sex stratified meta-analysis results

h2_liability <- h2_liability_varying_K_MA %>%
  filter((sex == "female" & K == 0.20) |
           (sex == "male" & K == 0.10)) %>%
  mutate(scale = "Liability") %>%
  select(-K)


h2_plot_sexes <- ggplot(h2_liability, aes(x = sex, y = mean_hsq * 100, colour = sex)) +
  geom_point() +
  geom_linerange(aes(ymin = hpdi_lower * 100, ymax = hpdi_upper * 100)) +
  geom_point(data = data.frame(sex = 'across', mean_hsq = NA)) +
  geom_linerange(data = data.frame(sex = 'across', dummy_x = 'male', ymin = -999, ymax = -999),
                 aes(x = dummy_x, ymin = ymin, ymax = ymax, colour = sex),
                 inherit.aes = FALSE, show.legend = TRUE) +
  scale_y_continuous(bquote('Liability ' ~{h^2} [SNP]~' (%)'),
                     limits = c(8, 12.5)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_manual(
    values = c('female' = '#FDE725FF', 
               'male' = '#440154FF', 
               'across' = '#20A387FF'),
    breaks = c("male", "female", "across"),
    labels = c('male' = "Male", 
               'female' = "Female", 
               'across' = "Across sex")
  ) +
  annotate("text", x = 1.5, y = 12.5, label = "P(F > M) = 100%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"))
h2_plot_sexes


### Polygenicity Plot ###

pi_plot <- ggplot(all_data_pi %>% filter(Type == "MA"), aes(x = sex, y = mean_pi, colour = sex)) +
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
  annotate("text", x = 1.5, y = 0.03, label = "P(F > M) = 100%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

pi_plot


### Selection Parameter Plot ###

S_plot <- ggplot(all_data_S %>% filter(Type == "MA"), aes(x = sex, y = mean_S, colour = sex)) +
  geom_point(position = position_dodge(width = -0.8)) +
  geom_linerange(aes(ymin = hpdi_lower, ymax = hpdi_upper),
                 position = position_dodge(width = -0.8)) +
  scale_x_discrete("",
                   limits = c("male", "female"),
                   labels = c("Male", "Female")) +
  scale_colour_viridis_d(end = 0, begin = 1) +
  guides(color = "none") +
  scale_y_continuous(expression("Selection parameter (" * italic(S) * ")"),
                     limits = c(-0.35, 0.2),
                     breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2),
                     labels = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2)) +
  annotate("text", x = 1.5, y = 0.19, label = "P(F > M) = 78%", size = 3) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

S_plot



### Female - male genetic correlation plot ###


rg_plot <- ggplot(all_data_rg, aes(x = estimate, y = test, colour = test)) +
  geom_vline(xintercept = 1,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point() +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) +
  scale_y_discrete("",
                   breaks = c("FF_across",
                              "MM_across",
                              "MF_across",
                              "MF_within",
                              "MF_MA"),
                   limits = c("FF_across",
                              "MM_across",
                              "MF_across",
                              "MF_within",
                              "MF_MA"),
                   labels = c("Female vs female \n(across cohorts)",
                              "Male vs male \n(across cohorts)",
                              "Male vs female \n(across cohorts)",
                              "Male vs female \n(within cohorts)",
                              "Male vs female \n(meta-analysed sumstats)")) +
  scale_colour_manual("Comparison",
                      values = c("FF_across" = "#FDE725FF",
                                 "MM_across" = "#440154FF",
                                 "MF_across" = "#20A387FF",
                                 "MF_within" = "#20A387FF",
                                 "MF_MA" = "#20A387FF"),
                      breaks = c("FF_across", "MM_across", "MF_across"),
                      labels = c("Female", "Male", "Across sex")) +
  guides(color = "none") +
  scale_x_continuous(expression(italic(r[g]))) +
  annotate("text", x = 0.97, y = "MF_MA", label = "*", size = 6) +
  annotate("text", x = 0.79, y = "MF_across", label = "*", size = 6) +
  annotate("text", x = 0.85, y = "FF_across", label = "*", size = 6) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

rg_plot



### Meta-Analysis of Female vs Male Betas Plot: Correlation (R) ###

R_plot <- ggplot(all_data_R, aes(x = estimate, y = test, colour = test)) +
  geom_vline(xintercept = 1,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point() +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) +
  scale_y_discrete("",
                   breaks = c("FF_r_across",
                              "MM_r_across",
                              "MF_r_across",
                              "MF_r_within",
                              "MF_MA"),
                   limits = c("FF_r_across",
                              "MM_r_across",
                              "MF_r_across",
                              "MF_r_within",
                              "MF_MA"),
                   labels = c("Female vs female (across cohorts)",
                              "Male vs male (across cohorts)",
                              "Male vs female (across cohorts)",
                              "Male vs female (within cohorts)",
                              "Male vs female (meta-analysed sumstats)")) +
  scale_colour_manual("Comparison",
                      breaks = c("FF_r_across",
                                 "MM_r_across",
                                 "MF_r_across",
                                 "MF_r_within",
                                 "MF_MA"),
                      values = c("#FDE725FF", "#440154FF", "#20A387FF", "#20A387FF", "#20A387FF")) +
  guides(color = "none") +
  scale_x_continuous("Correlation") +
  annotate("text", x = 0.87, y = "MF_MA", label = "*", size = 6) +
  annotate("text", x = 0.55, y = "MF_r_within", label = "*", size = 6) +
  annotate("text", x = 0.4, y = "MF_r_across", label = "*", size = 6) +
  annotate("text", x = 0.37, y = "MM_r_across", label = "*", size = 6) +
  annotate("text", x = 0.5, y = "FF_r_across", label = "*", size = 6) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

R_plot




### Mixer Venn ###

# Data: sizes of each group
female <- 6133
male <- 0
overlap <- 7111

sizes <- c(Female = female + overlap,
           Male = male + overlap)

# Calculate circle radii proportional to the square root of sizes
r_female <- sqrt(sizes["Female"] / pi)
r_male <- sqrt(sizes["Male"] / pi)

# Calculate offset to shift Male circle to the right so its right edge touches Female's right edge
offset <- r_female - r_male

# Create circle data with adjusted position for Male circle
circle_data <- data.frame(
  x = c(1, 1 + offset),  # Offset the Male circle to the right
  y = c(1, 1),           # Keep y-coordinates the same for female and male circle
  radius = c(r_female, r_male),
  group = c("Female", "Male")
)

# Plot the Euler diagram
mixer_venn <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6) +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  guides(fill = "none") +
  coord_fixed() +
  # labs(title = "Number of causal variants") +
  annotate("text", x = circle_data$x[1] - 50, y = circle_data$y[1], label = female, size = 4, color = "black", family = "Calibri") +
  annotate("text", x = circle_data$x[2], y = circle_data$y[2], label = overlap, size = 4, color = "black", family = "Calibri") +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

mixer_venn




### gwas-pw Venn ###

# Data: sizes of each group
female <- 3
male <- 0
overlap <- 42

sizes <- c(Female = female + overlap,
           Male = male + overlap)

# Calculate circle radii proportional to the square root of sizes
r_female <- sqrt(sizes["Female"] / pi)
r_male <- sqrt(sizes["Male"] / pi)

# Calculate offset to shift Male circle to the right so its right edge touches Female's right edge
offset <- r_female - r_male

# Create circle data with adjusted position for Male circle
circle_data <- data.frame(
  x = c(1, 1 + offset),  # Offset the Male circle to the right
  y = c(1, 1),           # Keep y-coordinates the same for female and male circle
  radius = c(r_female, r_male),
  group = c("Female", "Male")
)

# Plot the Euler diagram
gwaspw_venn <- ggplot(circle_data) +
  geom_circle(aes(x0 = x, y0 = y, r = radius, fill = group), alpha = 0.6) +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  guides(fill = "none") +
  coord_fixed() +
  # labs(title = "Number of genomic regions") +
  annotate("text", x = circle_data$x[1] - 4.5, y = circle_data$y[1], label = female, size = 4, color = "black", family = "Calibri") +
  annotate("text", x = circle_data$x[2], y = circle_data$y[2], label = overlap, size = 4, color = "black", family = "Calibri") +
  geom_segment(aes(x = circle_data$x[1] - 4.2, xend = circle_data$x[1] - 3.65,
                   y = circle_data$y[1], yend = circle_data$y[1]),
               color = "black", size = 0.5) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

gwaspw_venn


### Create Panelled Figure ###

figure <- (guide_area() /
  ((h2_plot_sexes | pi_plot | S_plot ) + plot_layout(widths = c(1, 1, 1))) /
  ((rg_plot | R_plot) + plot_layout(widths = c(1, 1))) /
  ((mixer_venn | gwaspw_venn) + plot_layout(widths = c(1, 1)))) +
  plot_layout(heights = c(0.2, 1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = "top")

figure

outfile <- paste(directory, "Compare_sexes_h2_pi_S_rg_R_mixer_gwaspw.png", sep="")
ggsave(figure, width = 21, height = 29, unit = "cm", file = outfile)

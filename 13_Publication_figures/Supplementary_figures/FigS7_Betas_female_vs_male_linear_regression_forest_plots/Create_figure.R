# Create publication figure of forest plots for male beta vs female beta linear regression slope and intercept

directory = "/path/10_Publication_figures/Supplementary_figures/Betas_female_vs_male_linear_regression_forest_plots/"



### Packages ###
library(patchwork)
library(scales)      # change scientific notation on plots
library(tidyverse)


### Load Data ###
# Linear regression data
load("/path/09_Effect_sizes_plots/Adams_SNPs/all_data_for_plot_for_slopes.RData")

load("/path/09_Effect_sizes_plots/Adams_SNPs/all_data_for_plot_for_intcpts.RData")


## Meta-Analysis of Female vs Male Betas Plot: Linear regression slope ###

Slope_plot <- ggplot(all_data_slope_plot, aes(x = mean_estimate, y = test, colour = test)) +
  geom_vline(xintercept = 1,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point() +
  geom_linerange(aes(xmin = mean_estimate_ci_lower, xmax = mean_estimate_ci_upper)) +
  scale_y_discrete("",
                   breaks = c("FF",
                              "MM",
                              "MF_slope_across",
                              "MF_slope_within",
                              "MF_MA"),
                   limits = c("FF",
                              "MM",
                              "MF_slope_across",
                              "MF_slope_within",
                              "MF_MA"),
                   labels = c("Female vs female \n(across cohorts)",
                              "Male vs male \n(across cohorts)",
                              "Male vs female \n(across cohorts)",
                              "Male vs female \n(within cohorts)",
                              "Male vs female \n(meta-analysed sumstats)")) +
  scale_colour_manual("Comparison",
                      breaks = c("FF",
                                 "MM",
                                 "MF_slope_across",
                                 "MF_slope_within",
                                 "MF_MA"),
                      values = c("#FDE725FF", "#440154FF", "#20A387FF", "#20A387FF", "#20A387FF")) +
  scale_x_continuous("Linear regression slope") +
  annotate("text", x = 0.9, y = "MF_MA", label = "*", size = 6) +
  annotate("text", x = 0.62, y = "MF_slope_within", label = "*", size = 6) +
  annotate("text", x = 0.54, y = "MF_slope_across", label = "*", size = 6) +
  annotate("text", x = 0.35, y = "MM", label = "*", size = 6) +
  annotate("text", x = 0.45, y = "FF", label = "*", size = 6) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

Slope_plot



### Meta-Analysis of Female vs Male Betas Plot: Linear regression intercept ###

Intcpt_plot <- ggplot(all_data_intcpt_plot, aes(x = mean_estimate, y = test, colour = test)) +
  geom_vline(xintercept = 0,
             colour = "grey30",
             linetype = 'dashed') +
  geom_point() +
  geom_linerange(aes(xmin = mean_estimate_ci_lower, xmax = mean_estimate_ci_upper)) +
  scale_y_discrete("",
                   breaks = c("FF",
                              "MM",
                              "MF_intcpt_across",
                              "MF_intcpt_within",
                              "MF_MA"),
                   limits = c("FF",
                              "MM",
                              "MF_intcpt_across",
                              "MF_intcpt_within",
                              "MF_MA"),
                   labels = c("Female vs female (across cohorts)",
                              "Male vs male (across cohorts)",
                              "Male vs female (across cohorts)",
                              "Male vs female (within cohorts)",
                              "Male vs female (meta-analysed sumstats)")) +
  scale_colour_manual("Comparison",
                      breaks = c("FF",
                                 "MM",
                                 "MF_intcpt_across",
                                 "MF_intcpt_within",
                                 "MF_MA"),
                      values = c("#FDE725FF", "#440154FF", "#20A387FF", "#20A387FF", "#20A387FF")) +
  scale_x_continuous("Linear regression intercept",
                     limits = c(-0.0006, 0.0015),
                     breaks = seq(-0.0006, 0.0015, 0.0004),
                     labels = function(x) ifelse(x == 0, "0", label_number(accuracy = 0.0001)(x))) +
  annotate("text", x = 0.001, y = "MF_intcpt_across", label = "*", size = 6) +
  annotate("text", x = 0.00065, y = "MM", label = "*", size = 6) +
  annotate("text", x = 0.00055, y = "FF", label = "*", size = 6) +
  theme_classic() +
  theme(text = element_text(family = "Calibri"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black"),
        legend.position = "none")

Intcpt_plot


### Create Panelled Figure ###

figure <- Slope_plot | Intcpt_plot+
  plot_annotation(tag_levels = 'A')


figure

outfile <- paste(directory, "Male_vs_female_betas_linear_regression.png", sep="")
ggsave(figure, width = 21, height = 12, unit = "cm", file = outfile)


args = commandArgs(trailingOnly=TRUE)
infile = args[1]
filename = args[2]
p_threshold = as.numeric(args[3])

#########################################################################

### Packages ###
library(data.table)   # To load data
library(qqman)        # To create QQ Plot
library(tidyverse)    

### Load Data ###
GWAS_df <- read.table(infile, header = TRUE, stringsAsFactors = F)

### Manhattan Plot ###
# Create cumulative BP position for x axis
cum_BP <- GWAS_df %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = dplyr::lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHR, bp_add)

Man_plot_df <- GWAS_df %>% 
  inner_join(cum_BP, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

# Find central BP for each chromosome for adding chromsome no. to x axis
axis_set <- Man_plot_df %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

# Set max y axis value
ylim <- Man_plot_df %>%
  filter(P == min(P)) %>%
  slice(1) %>%
  mutate(ylim = abs(floor(log10(P))) + 2) %>%
  pull(ylim)

# Set significance threshold
sig <- p_threshold

# Plot
Man_plot <- ggplot(Man_plot_df, 
       aes(x = bp_cum, y = -log10(P), color = as_factor(CHR), size = -log10(P))) +
  geom_point(size = 0.6, alpha = 0.75) +
  geom_hline(yintercept = -log10(sig), 
             color = "grey40",
             linetype = "dashed") +
  scale_x_continuous(label = axis_set$CHR, 
                     breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), 
                                  unique(length(axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = "Chromosome",
       y = expression(log[10](p-value))) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))

outfile <- paste(filename, ".Manhattan.png", sep="")
ggsave(Man_plot, width = 25, height = 10, unit = "cm", file = outfile)

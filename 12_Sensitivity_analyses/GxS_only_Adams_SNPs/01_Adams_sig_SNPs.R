# Make file containing Adams genome-wide significant SNPs
# Using Adams sumstats for all cohorts including 23&me, Europeans only

Adams <- read.table("/path/Summary_statistics_published/Depression_Adams_2024/including_23andMe/Adams_Full_EUR_2025_FromIMB.txt", header = T)

Adams_sig <- Adams %>% 
  filter(p < 5e-08) %>% 
  arrange(p)
# 42,317

write.table(Adams_sig, file = "/path/08_Sensitivity_analyses/Adams2024_topSNPs/Adams_Full_EUR_2025_FromIMB_sig_SNPs.txt", col.names = T, row.names = F, quote = F, sep = "\t")

## load power simulation output to compute power for one sample and effect size

# required packages
library(tidyverse)

# column types of power simulation files
sim_col_types <- cols(
  iteration = col_double(),
  perturbation = col_character(),
  perturbed = col_logical(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  disp_outlier_deseq2 = col_logical()
)

# load power simulation output files
power_sim <- lapply(snakemake@input, FUN = read_csv, col_types = sim_col_types, progress = FALSE)

# combine into one data frame
power_sim <- bind_rows(power_sim)

# compute power to detect repression effect of each enhancer - gene pair
power <- power_sim %>% 
  filter(perturbed == TRUE, !is.na(logFC)) %>% 
  group_by(iteration) %>% 
  mutate(pval_adj = p.adjust(pvalue, method = "fdr")) %>% 
  group_by(perturbation, gene, disp_outlier_deseq2) %>%
  summarize(power = mean(pval_adj < snakemake@params$fdr & logFC < 0), .groups = "drop") %>% 
  arrange(desc(power))

# add effect simulated size
power$effect_size <- as.numeric(snakemake@wildcards$effect)

# save power calculations to output file
write_csv(power, file = snakemake@output[[1]])

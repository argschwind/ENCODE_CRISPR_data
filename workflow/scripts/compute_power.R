## load power simulation output to compute power for one sample and effect size

# save.image("compute_pwr.rda")
# stop()

library(tidyverse)

# get vector of power simulation output files and set names to iteration number
power_sim_files <- unlist(snakemake@input)
names(power_sim_files) <- sub(".+_rep([[:digit:]]+)_.+", "\\1", power_sim_files)

# column types of power simulation files
col_types <- cols(
  iteration = col_double(),
  perturbation = col_character(),
  pert_gene = col_character(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  disp_outlier_deseq2 = col_logical()
)

# load power simulation output files
power_sim <- lapply(power_sim_files, FUN = read_csv, col_types = col_types, progress = FALSE)

# combine into one data frame
power_sim <- bind_rows(power_sim)

# compute power to detect repression effect of each enhancer - gene pair
#power_pert_gene <- power_sim %>% 
#  group_by(iteration, pert_gene) %>% 
#  mutate(pval_adj = p.adjust(pvalue, method = "fdr")) %>% 
#  filter(pert_gene == gene) %>% 
#  group_by(perturbation, gene) %>% 
#  summarize(power = mean(pval_adj < 0.05)) %>% 
#  arrange(desc(power))

# compute power to detect repression effect of each enhancer - gene pair
power <- power_sim %>% 
  filter(pert_gene == gene) %>% 
  group_by(iteration) %>% 
  mutate(pval_adj = p.adjust(pvalue, method = "fdr")) %>% 
  group_by(perturbation, gene, disp_outlier_deseq2) %>%
  summarize(power = mean(pval_adj < snakemake@params$fdr), .groups = "drop") %>% 
  arrange(desc(power))

# save power calculations to output file
write_csv(power, path = snakemake@output[[1]])

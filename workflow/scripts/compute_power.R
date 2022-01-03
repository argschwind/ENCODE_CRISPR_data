## load power simulation output to compute power for one sample and effect size

library(tidyverse)

# column types of power simulation files
col_types <- cols(
  perturbation = col_character(),
  iteration = col_double(),
  pert_gene = col_character(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  disp_outlier_deseq2 = col_logical()
)

# load power simulation output
power_sim <- read_csv(snakemake@input[[1]], col_types = col_types)

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
  group_by(perturbation, gene) %>% 
  summarize(power = mean(pval_adj < snakemake@params$fdr)) %>% 
  arrange(desc(power))

# add DESeq2 dispersion outlier information
power <- power_sim %>% 
  select(gene, disp_outlier_deseq2) %>% 
  distinct() %>% 
  left_join(x = power, y = ., by = "gene")

# save power calculations to output file
write_csv(power, path = snakemake@output[[1]])

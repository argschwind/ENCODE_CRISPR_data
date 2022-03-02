## Combine differential expression and power simulation results to create output data sets

# required packages
library(dplyr)
library(tidyr)
library(readr)

# load input data ----------------------------------------------------------------------------------

# column types in DE results 
diff_expr_cols <- cols(
  perturbation = col_character(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  pval_adj = col_double(),
  cells = col_double(),
  avg_expr = col_double(),
  pert_level = col_character(),
  chr = col_character(),
  pert_start = col_double(),
  pert_end = col_double(),
  gene_tss = col_double(),
  gene_strand = col_character(),
  dist_to_tss = col_double()
)

# load differential expression results
diff_expr <- read_tsv(snakemake@input$diff_expr, col_types = diff_expr_cols, progress = FALSE)

# column types in power files
power_cols <- cols(
  perturbation = col_character(),
  gene = col_character(),
  disp_outlier_deseq2 = col_logical(),
  power = col_double(),
  effect_size = col_double()
)

# load power files and combine into one data frame
power <- snakemake@input$power_sim %>% 
  lapply(FUN = read_tsv, col_types = power_cols, progress = FALSE) %>% 
  bind_rows()

# add simulated power to differential expression results -------------------------------------------

# reformat effect size to percentage and transform to wide format with new column names
power <- power %>%
  mutate(effect_size = effect_size * 100) %>% 
  pivot_wider(names_from = effect_size, names_prefix = "power_effect_size_", values_from = power)

# add power to DE results
output <- left_join(diff_expr, power, by = c("perturbation", "gene"))

# write to .tsv file
write_tsv(output, file = snakemake@output[[1]])

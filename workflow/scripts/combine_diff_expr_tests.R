## Combine per-chromosome differential expression tests into one file

# required packages
library(dplyr)
library(readr)

# column types in DE result files
de_cols <- cols(
  perturbation = col_character(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  pval_adj = col_double(),
  pert_chr = col_character(),
  pert_start = col_integer(),
  pert_end = col_integer(),
  gene_chr = col_character(),
  gene_tss = col_integer(),
  gene_strand = col_character(),
  dist_to_tss = col_double(),
  pert_level = col_character(),
  target_type = col_character(),
  cells = col_double(),
  avg_expr = col_double()
)

# load DE results for each chromosome
de_results <- lapply(unlist(snakemake@input), FUN = read_tsv, col_types = de_cols)

# combine into one table and re-calculate multiple testing correction
output <- de_results %>% 
  bind_rows() %>% 
  mutate(pval_adj = p.adjust(pvalue, method = snakemake@params$padj_method))

# save to output file
write_tsv(output, file = snakemake@output[[1]])

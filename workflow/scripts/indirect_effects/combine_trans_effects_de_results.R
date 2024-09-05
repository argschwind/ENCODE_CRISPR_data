## Combine results from trans-acting differential expression test across perturbation batches and
## assign significance based on 5% FDR threshold in cis analysis

# save.image("combine.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# load all trans-effect results and combine into one table
results <- snakemake@input$trans_results %>% 
  lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = "batch")

# load cis differential expression output
cis_results <- read_tsv(snakemake@input$cis_results, show_col_types = FALSE)

# extract nominal p-value corresponding to 5% FDR on cis tests
pval_threshold_cis <- cis_results %>% 
  filter(pval_adj < 0.05) %>% 
  pull(pvalue) %>% 
  max()

# add significance to trans-effect results based on cis significance threshold
results <- results %>% 
  mutate(significant = pvalue <= pval_threshold_cis,
         regulated = significant & logFC < 0)

# write output to file
write_tsv(results, snakemake@output[[1]])

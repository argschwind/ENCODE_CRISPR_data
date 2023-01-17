## Combine main output files for both TAP-seq experiments into one file. This also filters out any
## duplicated control enhancer perturbations

# save.image("combine_de.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# column in results files
results_cols <- cols(
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

# load results for TAP-seq experiments
results_files <- c(chr8 = snakemake@input$chr8, chr11 = snakemake@input$chr11)
results <- lapply(results_files, FUN = read_tsv, col_types = results_cols, progress = FALSE)

# combine into one data frame
results <- bind_rows(results, .id = "sample")

# re-compute FDR
results <- mutate(results, pval_adj = p.adjust(pvalue, method = snakemake@params$p_adj_method))

# remove out of region controls for MYC, ZFPM2 and HB* and only keep chr11 control for GATA1
results <- results %>% 
  filter( (pert_chr == sample & gene_chr == sample) |
            (gene == "GATA1" & grepl(perturbation, pattern = "GATA1") & sample == "chr11") )

# remove sample column and sort output
results <- results %>% 
  select(-sample) %>% 
  arrange(pert_chr, pert_start, pert_end, gene_chr, gene_tss)

# write merged data to output file
write_tsv(results, file = snakemake@output[[1]])

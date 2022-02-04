## Estimate guide-guide variability from per element and per guide differential expression results.

# required packages and functions
library(tidyverse)
library(cowplot)
source(file.path(snakemake@scriptdir, "R_functions/estimate_guide_variability_fun.R"))

# column types in differential expression results
de_results_cols <- cols(
  perturbation = col_character(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  pval_adj = col_double(),
  cells = col_double(),
  avg_expr = col_double(),
  chr = col_character(),
  pert_start = col_double(),
  pert_end = col_double(),
  gene_tss = col_double(),
  gene_strand = col_character(),
  dist_to_tss = col_double()
)

# load differential expression results
per_target <- read_csv(snakemake@input$per_target, col_types = de_results_cols)
per_guide  <- read_csv(snakemake@input$per_guide,  col_types = de_results_cols)

# column types in guide targets
guide_targets_cols <- cols(
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  name = col_character(),
  strand = col_character(),
  spacer = col_character(),
  target_chr = col_character(),
  target_start = col_double(),
  target_end = col_double(),
  target_name = col_character(),
  target_strand = col_character()
)

# load guide target annotations
guide_targets <- read_tsv(snakemake@input$guide_targets, col_types = guide_targets_cols)

# rename column containing effect size in differential expression results
colnames(per_target)[colnames(per_target) == snakemake@params$effect_size_col] <- "effect_size"
colnames(per_guide)[colnames(per_guide)  == snakemake@params$effect_size_col] <- "effect_size"

# extract significant pairs from on per-target results and filter out pairs that increase expression
sig_per_target <- filter(per_target, pval_adj < snakemake@params$fdr_sig)
sig_per_target <- filter(sig_per_target, effect_size < 0)

# estimate guide variability from significant pairs
guide_var <- estimate_guide_variability(sig_per_target,
                                        per_guide = per_guide,
                                        guide_targets = guide_targets,
                                        cells_per_guide = snakemake@params$cells_per_guide,
                                        effect_size = snakemake@params$effect_size_col,
                                        pct_change_range = c(-1, 0),
                                        return_plots = TRUE)

# save guide variability and distribution fit to output files
write_csv(guide_var$guide_variability, file = snakemake@output$guide_var)
write_csv(guide_var$distribution, file = snakemake@output$distr_fit)

# arrange plots into one panel and save to file
plots <- plot_grid(plotlist = guide_var$plots, rel_widths = c(4, 3))
ggsave(plots, file = snakemake@output$plots, width = 8, height = 4)

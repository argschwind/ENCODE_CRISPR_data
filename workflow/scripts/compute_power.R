## Load power simulation output to compute power for one sample and effect size

# save.image("compute_pwr.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

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

# load power simulation output files, combine into one data frame and only retain results for
# simulated perturbations
power_sim <- snakemake@input$sim %>% 
  lapply(FUN = read_tsv, col_types = sim_col_types, progress = FALSE) %>% 
  bind_rows() %>% 
  filter(perturbed == TRUE, !is.na(logFC)) 

# column types in differential expression results file
real_col_types <- cols(
  perturbation = col_character(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  pval_adj = col_double(),
  pert_chr = col_character(),
  pert_start = col_double(),
  pert_end = col_double(),
  gene_chr = col_character(),
  gene_tss = col_double(),
  gene_strand = col_character(),
  dist_to_tss = col_double(),
  pert_level = col_character(),
  target_type = col_character(),
  cells = col_double(),
  avg_expr = col_double()
)

# load differential expression testing results
diff_expr <- read_tsv(snakemake@input$real, col_types = real_col_types, progress = FALSE)

# compute power to detect repression effect of each enhancer - gene pair depending on power
# calculation strategy
if (snakemake@params$multi_test_correction == "sim") {
  
  # compute power using FDR correction applied in simulated results
  power <- power_sim %>% 
    group_by(iteration) %>% 
    mutate(pval_adj = p.adjust(pvalue, method = snakemake@params$p_adj_method)) %>% 
    group_by(perturbation, gene, disp_outlier_deseq2) %>%
    summarize(power = mean(pval_adj < snakemake@params$pval_threshold & logFC < 0),
              .groups = "drop") %>% 
    arrange(desc(power))
  
} else if (snakemake@params$multi_test_correction == "real") {
  
  # get nominal p-value threshold for FDR corrected hits from real differential expression results
  pval_threshold <- max(pull(filter(diff_expr, pval_adj < snakemake@params$pval_threshold), pvalue))
  
  # compute power using nominal p-value threshold from real data
  power <- power_sim %>% 
    group_by(perturbation, gene, disp_outlier_deseq2) %>%
    summarize(power = mean(pvalue < pval_threshold & logFC < 0), .groups = "drop") %>% 
    arrange(desc(power))
  
} else {
  stop("Incorrect 'multi_test_correction' argument. Needs to be one of 'real' or 'sim'",
       call. = FALSE)
}

# add simulated effect size
power$effect_size <- as.numeric(snakemake@wildcards$effect)

# save power calculations to output file
write_tsv(power, file = snakemake@output[[1]])

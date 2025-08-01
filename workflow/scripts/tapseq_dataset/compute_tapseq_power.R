## Load power simulation output to compute power for both TAP-seq experiments. This also filters out
## any duplicated control enhancer perturbations

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

# get input files for all chromosomes
chr_input_files <- snakemake@input[grepl(names(snakemake@input), pattern = "^chr")]

# load power simulation output files for all chromosomes
power_sim <- lapply(chr_input_files, FUN = function (chr) {
  lapply(chr, FUN = read_tsv, col_types = sim_col_types, progress = FALSE)
})

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

# combine simulation results into one data frame and only retain results for simulated perturbations
power_sim <- power_sim %>% 
  lapply(bind_rows) %>% 
  bind_rows(.id = "sample") %>% 
  filter(perturbed == TRUE, !is.na(logFC))

# compute power to detect repression effect of each enhancer - gene pair depending on power
# calculation strategy
if (snakemake@params$multi_test_correction == "sim") {
  
  # compute power using FDR correction applied in simulated results
  power <- power_sim %>% 
    group_by(sample, iteration) %>% 
    mutate(pval_adj = p.adjust(pvalue, method = snakemake@params$p_adj_method)) %>% 
    group_by(sample, perturbation, gene, disp_outlier_deseq2) %>%
    summarize(power = mean(pval_adj < snakemake@params$pval_threshold & logFC < 0),
              .groups = "drop") %>% 
    arrange(sample, desc(power))
  
} else if (snakemake@params$multi_test_correction == "real") {
  
  # get nominal p-value threshold for FDR corrected hits from real differential expression results
  pval_threshold <- max(pull(filter(diff_expr, pval_adj < snakemake@params$pval_threshold), pvalue))
  
  # compute power using nominal p-value threshold from real data
  power <- power_sim %>% 
    group_by(sample, perturbation, gene, disp_outlier_deseq2) %>%
    summarize(power = mean(pvalue < pval_threshold & logFC < 0), .groups = "drop") %>% 
    arrange(sample, desc(power))
  
} else {
  stop("Incorrect 'multi_test_correction' argument. Needs to be one of 'real' or 'sim'",
       call. = FALSE)
}

# add column specifying whether a pair is an out-of-region control (include GATA1 in chr11)
power <- power %>% 
  mutate(oor_ctrl = case_when(
    sample == "chr8" & perturbation %in% c("HS2", "GATA1") ~ TRUE,
    sample == "chr8" & gene %in% c("HBD", "HBB", "HBE1", "GATA1") ~ TRUE,
    sample == "chr11" & perturbation %in% c("MYC", "ZFPM2") ~ TRUE,
    sample == "chr11" & gene %in% c("MYC", "ZFPM2") ~ TRUE,
    TRUE ~ FALSE
  ))

# remove pairs labelled as out-of-regions controls and remove added columns
power <- power %>% 
  filter(oor_ctrl == FALSE) %>% 
  select(-c(sample, oor_ctrl))

# add simulated effect size
power$effect_size <- as.numeric(snakemake@wildcards$effect)

# save power calculations to output file
write_tsv(power, file = snakemake@output[[1]])

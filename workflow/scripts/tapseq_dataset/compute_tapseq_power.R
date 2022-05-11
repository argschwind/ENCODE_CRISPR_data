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

# combine into one data frame
power_sim <- power_sim %>% 
  lapply(bind_rows) %>% 
  bind_rows(.id = "sample")

# compute power to detect repression effect of each enhancer - gene pair
power <- power_sim %>% 
  filter(perturbed == TRUE, !is.na(logFC)) %>% 
  group_by(sample, iteration) %>% 
  mutate(pval_adj = p.adjust(pvalue, method = snakemake@params$p_adj_method)) %>% 
  group_by(sample, perturbation, gene, disp_outlier_deseq2) %>%
  summarize(power = mean(pval_adj < snakemake@params$pval_threshold & logFC < 0),
            .groups = "drop") %>% 
  arrange(sample, desc(power))

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

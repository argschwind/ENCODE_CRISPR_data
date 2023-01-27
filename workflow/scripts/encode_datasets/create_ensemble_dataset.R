## Create combined CRISPR dataset for distal regulation CRISPR E-P benchmarking

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

## Load CRISPR data --------------------------------------------------------------------------------

# get all input files
infiles <- snakemake@input[names(snakemake@input) != ""]

# column types in ENCODE CRISPR data
crispr_data_cols <- cols(
  .default = col_double(),  # for power columns, which can vary in names
  chrom = col_character(),
  chromStart = col_integer(),
  chromEnd = col_integer(),
  name = col_character(),
  EffectSize = col_double(),
  strandPerturbationTarget = col_character(),
  PerturbationTargetID = col_character(),
  chrTSS = col_character(),
  startTSS = col_integer(),
  endTSS = col_integer(),
  strandGene = col_character(),
  EffectSize95ConfidenceIntervalLow = col_double(),
  EffectSize95ConfidenceIntervalHigh = col_double(),
  measuredGeneSymbol = col_character(),
  measuredEnsemblID = col_character(),
  guideSpacerSeq = col_character(),
  guideSeq = col_character(),
  Significant = col_logical(),
  pValue = col_double(),
  pValueAdjusted = col_double(),
  ValidConnection = col_character(),
  Notes = col_character(),
  Reference = col_character()
)

# load all CRISPR datasets
crispr_data <- lapply(infiles, FUN = read_tsv, col_types = crispr_data_cols, progress = FALSE)

## Pre-process input data --------------------------------------------------------------------------

# helper function to convert log-fold changes to percentage changes
logFC_to_pctChange <- function(logFC, base = exp(1)) {
  pct_change <- base ^ logFC - 1
  return(pct_change)
}

# function to convert all effect size columns in dataset to percent changes
covert_all_to_pctChange <- function(dat, effect_size) {
  if (effect_size != "pctChange") {
    base <- switch(effect_size, "logFC" = exp(1), "log2FC" = 2)
    dat <- dat %>% 
      mutate(EffectSize = logFC_to_pctChange(EffectSize, base = base),
             EffectSize95ConfidenceIntervalLow = logFC_to_pctChange(EffectSize95ConfidenceIntervalLow, base = base),
             EffectSize95ConfidenceIntervalHigh = logFC_to_pctChange(EffectSize95ConfidenceIntervalHigh, base = base))
  }
  return(dat)
}

# convert effect sizes to percent change in expression for datasets with log FC effect sizes
effect_sizes <- unlist(snakemake@params$effect_size)
crispr_data <- mapply(dat = crispr_data, effect_size = effect_sizes[names(crispr_data)],
                      FUN = covert_all_to_pctChange, SIMPLIFY = FALSE)

# combine into one data frame
crispr_data <- bind_rows(crispr_data, .id = "Dataset")

# filter all datasets for valid connections
crispr_data <- crispr_data %>% 
  filter(ValidConnection == "TRUE") %>% 
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)

## Combine CRISPR datasets -------------------------------------------------------------------------

# add unique identifier to each pair
crispr_data <- mutate(crispr_data,
                      pair_uid = paste(Dataset, name, sep = "|"))

# create GRanges object containing all enhancer-gene pairs for each dataset
eg_pairs <- with(crispr_data,
                 GRanges(seqnames = paste0(chrom,":", measuredGeneSymbol),
                         ranges = IRanges(chromStart, chromEnd), pair_uid = pair_uid))

# get unique enhancer gene pairs
merged_pairs <- reduce(eg_pairs, with.revmap = TRUE)

# add unique identifier for merged enhancer-gene pairs
merged_pairs$merged_uid <- seq_along(merged_pairs)

# convert to data frame and get pairs that need to be merged for each merged pair
merged_pairs <- merged_pairs %>% 
  as.data.frame() %>% 
  unnest(cols = revmap)

# add pair uid
merged_pairs$pair_uid <- eg_pairs$pair_uid[merged_pairs$revmap]

# add merged uid and merged coordinates to crispr_data
crispr_data <- merged_pairs %>% 
  select(merged_uid, pair_uid, merged_start = start, merged_end = end) %>% 
  left_join(x = crispr_data, y = ., by = "pair_uid")

# function to process overlapping pairs
process_overlapping_pairs <- function(pair) {
  
  # only retain significant pairs if there are any significants
  if (any(pair$Significant == TRUE)) {
    pair <- filter(pair, Significant == TRUE)
  }
  
  # pick top powered pair if more than one remain (e.g. when all overlapping pairs were negative)
  slice_max(pair, n = 1, order_by = PowerAtEffectSize25, with_ties = FALSE)
  
}

# split crispr_data by merged_uid
crispr_data_split <- split(crispr_data, f = crispr_data$merged_uid)

# process overlapping pairs
output <- lapply(crispr_data_split, FUN = process_overlapping_pairs)
output <- bind_rows(output)

# Create output files ------------------------------------------------------------------------------

# rearrange columns for output in ENCODE format
encode <- output %>% 
  select(all_of(names(crispr_data_cols$cols)), starts_with("PowerAtEffectSize"), Dataset) %>% 
  relocate(starts_with("PowerAtEffectSize"), .before = ValidConnection)

# sort according to genomic coordinates
encode <- arrange(encode, chrom, chromStart, chromEnd, measuredGeneSymbol)

# write to output files
write_tsv(encode, file = snakemake@output[[1]])

## Filter HCR-FlowFISH data

# save.image("filt.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# column types in ENCODE CRISPR data
crispr_data_cols <- cols(
  .default = col_guess(),  # for optional columns, which can vary in names and format
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
  Reference = col_character(),
  distToTSS = col_double()
)

# load CRISPR dataset
dat <- read_tsv(snakemake@input[[1]], col_types = crispr_data_cols, progress = FALSE)

# filter based on distance to TSS
dist_filter <- as.integer(snakemake@params$tss_to_dist)
dat <- filter(dat, abs(distToTSS) > dist_filter[[1]], abs(distToTSS) <= dist_filter[[2]])

# write to output
write_tsv(dat, file = snakemake@output[[1]])

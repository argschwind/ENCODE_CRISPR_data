## Extract columns for EP benchmarking CRISPR data format from full ENCODE data and filter for
## valid connections, distance to TSS and power.

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# helper function to convert log-fold changes to percentage changes
logFC_to_pctChange <- function(logFC, base = exp(1)) {
  pct_change <- base ^ logFC - 1
  return(pct_change)
}

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
  Reference = col_character()
)

# load CRISPR dataset
dat <- read_tsv(snakemake@input[[1]], col_types = crispr_data_cols, progress = FALSE)

# parse effect_size argument
effect_size <- match.arg(snakemake@params$effect_size, choices = c("pctChange", "logFC", "log2FC"))

# convert effect size to 'percent change in expression' if required
if (effect_size != "pctChange") {
  message("Converting EffectSize from ", effect_size, " to pctChange.")
  base <- switch(effect_size, "logFC" = exp(1), "log2FC" = 2)
  dat <- mutate(dat, EffectSize = logFC_to_pctChange(EffectSize, base = base))
}

# add 'Regulated' column, labeling detected enhancer-gene pairs
if (is.null(snakemake@params$min_pct_change)) {
  dat$Regulated <- dat$Significant == TRUE & dat$EffectSize < 0
} else {
  dat$Regulated <- dat$Significant == TRUE & dat$EffectSize <= snakemake@params$min_pct_change
}

# add cell type column
dat <- mutate(dat, CellType = snakemake@params$cell_type)

# extract required columns for EP benchmarking data and sort according to coordinates
output <- dat %>%
  select(chrom, chromStart, chromEnd, name, EffectSize, chrTSS, startTSS,
         endTSS, measuredGeneSymbol, Significant, pValueAdjusted, starts_with("PowerAtEffectSize"),
         ValidConnection, CellType, Reference, Regulated, any_of("Dataset")) %>% 
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)

# write output to file
write_tsv(output, file = snakemake@output[[1]])

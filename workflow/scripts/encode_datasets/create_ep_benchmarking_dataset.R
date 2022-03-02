## Extract columns for EP benchmarking CRISPR data format from full ENCODE data and filter for
## valid connections, distance to TSS and power.

# save.image("ep_input.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

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
  guideSpacerSeq = col_logical(),
  guideSeq = col_character(),
  Significant = col_logical(),
  pValue = col_double(),
  pValueAdjusted = col_double(),
  ValidConnection = col_character(),
  CellType = col_character(),
  Reference = col_character(),
  distToTSS = col_double(),
  avgGeneExpr = col_double(),
  nPertCells = col_double()
)

# load CRISPR dataset
dat <- read_tsv(snakemake@input[[1]], col_types = crispr_data_cols, progress = FALSE)

# distance to TSS filter (min and max)
dist_filter <- as.integer(snakemake@params$tss_to_dist)

# filter for distance to TSS
dat <- filter(dat, abs(distToTSS) > dist_filter[[1]], abs(distToTSS) <= dist_filter[[2]])

# extract required columns for EP benchmarking data and sort according to coordinates
output <- dat %>%
  select(chrom, chromStart, chromEnd, name, EffectSize, chrTSS, startTSS,
         endTSS, measuredGeneSymbol, Significant, pValueAdjusted, starts_with("PowerAtEffectSize"),
         ValidConnection, CellType, Reference) %>% 
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)

# write output to file
write_tsv(output, file = snakemake@output[[1]])

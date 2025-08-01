## Split Nasser et al., 2021 data into different benchmarking datasets based on cell types

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# column types in Nasser2021 crispr dataset
nasser_data_cols <- cols(
  chrom = col_character(),
  chromStart = col_double(),
  chromEnd = col_double(),
  name = col_character(),
  EffectSize = col_double(),
  chrTSS = col_character(),
  startTSS = col_double(),
  endTSS = col_double(),
  measuredGeneSymbol = col_character(),
  Significant = col_logical(),
  pValueAdjusted = col_double(),
  PowerAtEffectSize25 = col_double(),
  ValidConnection = col_character(),
  CellType = col_character(),
  Reference = col_character(),
  Regulated = col_logical()
)

# load Nasser2021 crispr dataset
dat <- read_tsv(snakemake@input[[1]], col_types = nasser_data_cols, progress = FALSE)

# subset data to selected cell types
dat <- filter(dat, CellType %in% snakemake@params$subset)

# write subset to output files
write_tsv(dat, file = snakemake@output[[1]])

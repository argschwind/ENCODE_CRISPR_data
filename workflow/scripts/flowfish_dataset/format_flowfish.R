## Reformat Engreitz lab FlowFISH data into EP benchmarking CRISPR data file format

# required packages
library(dplyr)
library(readr)

# column types input flowfish data
flowfish_cols <- cols(
  chrPerturbationTarget = col_character(),
  startPerturbationTarget = col_double(),
  endPerturbationTarget = col_double(),
  chrTSS = col_character(),
  startTSS = col_double(),
  endTSS = col_double(),
  GeneSymbol = col_character(),
  CellType = col_character(),
  name = col_character(),
  Reference = col_character(),
  EffectSize = col_double(),
  Significant = col_logical(),
  padj = col_double(),
  Regulated = col_logical(),
  IncludeInModel = col_logical(),
  Reason = col_character(),
  expt.class = col_character(),
  PowerAtEffectSize.0.25 = col_double(),
  PerturbMethod = col_character(),
  RNAReadoutMethod = col_character()
)

# load input data
dat <- read_tsv(snakemake@input[[1]], col_types = flowfish_cols, progress = FALSE)

# rename power column(s)
colnames(dat) <- sub("PowerAtEffectSize.0.", "PowerAtEffectSize", colnames(dat))

# create ValidConnection column
dat <- dat %>% 
  mutate(ValidConnection = if_else(IncludeInModel == TRUE, true = "TRUE", false = Reason))

# create name column based on enhancer coordinates and target gene
dat <- mutate(dat, name = paste0(GeneSymbol, "|", chrPerturbationTarget, ":",
                                   startPerturbationTarget, "-", endPerturbationTarget, ":*"))

# extract required columns for EP benchmarking format and sort according to enhancer coordinates
output <- dat %>% 
  select(chrom = chrPerturbationTarget, chromStart = startPerturbationTarget,
         chromEnd = endPerturbationTarget, name, EffectSize, chrTSS, startTSS, endTSS,
         measuredGeneSymbol = GeneSymbol, Significant, pValueAdjusted = padj,
         PowerAtEffectSize25, ValidConnection, CellType, Reference, Regulated) %>% 
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)

# write output to file
write_tsv(output, file = snakemake@output[[1]])

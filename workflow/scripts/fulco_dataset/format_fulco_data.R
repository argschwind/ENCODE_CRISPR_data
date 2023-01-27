## Reformat Engreitz lab Fulco et al., 2019 data into ENCODE CRISPR data file format

# required packages
library(dplyr)
library(readr)

## Load input data ---------------------------------------------------------------------------------

# column types input Fulco data
fulco_cols <- cols(
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

# load FlowFish input data
dat <- read_tsv(snakemake@input$data, col_types = fulco_cols, progress = FALSE)

# column types input tss annotations
tss_cols <- cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character()
)

# load TSS annotations
tss <- read_tsv(snakemake@input$tss, col_types = tss_cols, col_names = names(tss_cols$cols))

## Reformat and add columns ------------------------------------------------------------------------

# rename power column(s)
colnames(dat) <- sub("PowerAtEffectSize.0.", "PowerAtEffectSize", colnames(dat))

# create ValidConnection column
dat <- dat %>% 
  mutate(ValidConnection = if_else(IncludeInModel == TRUE, true = "TRUE", false = Reason))

# add missing perturbation and name columns
dat <- dat %>% 
  mutate(strandPerturbationTarget = ".",
         PerturbationTargetID = paste0(chrPerturbationTarget, ":", startPerturbationTarget, "-",
                                       endPerturbationTarget, ":", strandPerturbationTarget),
         name = paste0(GeneSymbol, "|", PerturbationTargetID))

# add target gene strand information
dat <- dat %>% 
  mutate(gene = sub("-TSS1$", "", GeneSymbol)) %>%  # strip TSS information from some genes symbols
  left_join(select(tss, gene = name, strandGene = strand), by = "gene") %>% 
  select(-gene)

# reformat references to prettier names
dat <- mutate(dat, Reference = sub("([[:alpha:]]+)([[:digit:]]+)", "\\1 et al., \\2", Reference))

# add miscellaneous columns
dat <- dat %>% 
  mutate(EffectSize95ConfidenceIntervalLow = NA_real_,
         EffectSize95ConfidenceIntervalHigh = NA_real_,
         measuredEnsemblID = NA_character_, guideSpacerSeq = NA_character_,
         guideSeq = NA_character_, pValue = NA_real_, Notes = NA_character_)

# extract columns for ENCODE CRISPR element quantification format
output <- dat %>% 
  select(chrom = chrPerturbationTarget, chromStart = startPerturbationTarget,
         chromEnd = endPerturbationTarget, name, EffectSize, strandPerturbationTarget, 
         PerturbationTargetID, chrTSS, startTSS, endTSS, strandGene,
         EffectSize95ConfidenceIntervalLow, EffectSize95ConfidenceIntervalHigh,
         measuredGeneSymbol = GeneSymbol, measuredEnsemblID, guideSpacerSeq, guideSeq,
         Significant, pValue, pValueAdjusted = padj, starts_with("PowerAtEffectSize"),
         ValidConnection, Notes, Reference)

# write output to file
write_tsv(output, file = snakemake@output[[1]])

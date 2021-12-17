## Create bed-like file containing re-mapped Gasperini guide positions

# required packages
library(tidyverse)

# column types in remapped guides file
guides_cols <- cols(
  GuideID = col_character(),
  Spacer = col_character(),
  Target_Site = col_character(),
  chr.candidate_enhancer = col_character(),
  start.candidate_enhancer = col_integer(),
  stop.candidate_enhancer = col_integer(),
  Category = col_character(),
  original.order = col_integer(),
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  score = col_double(),
  strand = col_character(),
  nPerfectMatches = col_integer()
)

# load re-mapped guides
guides <- read_tsv(snakemake@input$guides, col_types = guides_cols)

# load chromosome sizes file to get order for sorting (like bedtools sort)
chrs_cols <- cols(chr = col_character(), length = col_integer())
chrs <- read_tsv(snakemake@input$chrs, col_types = chrs_cols, col_names = names(chrs_cols$cols))

# filter out non-targeting controls and extract guide position columns in bed-like format
guides_bed <- guides %>% 
  filter(!is.na(start), Category != "NTC") %>% 
  select(chr, start, end, GuideID, score, strand, Spacer)

# add G to 19 bp long spacer sequences
guides_bed <- guides_bed %>% 
  mutate(Spacer = if_else(nchar(Spacer) == 19, true = paste0(Spacer, "G"), false = Spacer))

# sort according to chromosome (from chrs file) and position
guides_bed <- guides_bed %>% 
  mutate(chr = factor(chr, levels = chrs$chr)) %>% 
  arrange(chr, start, end, GuideID)

# write to output file
write_tsv(guides_bed, file = snakemake@output[[1]], col_names = FALSE)
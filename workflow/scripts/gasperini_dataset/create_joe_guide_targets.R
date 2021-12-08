## Create a guide targets file from Joe's approach of merging proximal re-mapped guides into regions

# required packages
library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)

# columns types in the re-mapped gRNAs bed file
guide_cols <- cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character(),
  spacer = col_character()
)

# read file containing re-mapped guides and their target regions
guides <- read_tsv(snakemake@input[[1]], col_types = guide_cols, col_names = names(guide_cols$cols))

# create GRanges from re-mapped guide coordinates
guides_gr <- makeGRangesFromDataFrame(guides, keep.extra.columns = TRUE,
                                      starts.in.df.are.0based = TRUE)

# merge nearby guides to create target regions
targets <- reduce(guides_gr + snakemake@params$extend_bp, with.revmap = TRUE, ignore.strand = TRUE)

# create tibble with target coordinates in bed format (convert start coordinates back to 0-based)
# as.data.frame conversion does not work with revmap list
targets_tbl <- tibble(
    target_chr = as.character(seqnames(targets)), target_start = start(targets),
    target_end = end(targets), target_strand = as.character(strand(targets)),
    revmap = as.list(targets$revmap) ) %>% 
  mutate(target_start = target_start - 1,
         target_name = paste0(target_chr, ":", target_start, "-", target_end),
         target_score = 666) %>% 
  select(target_chr, target_start, target_end, target_name, target_score, target_strand, revmap)

# flatten revmap column into a regular column using unnest
targets_tbl <- unnest(targets_tbl, cols = revmap)

# merge with guides based on revmap column (index of guides merged for a given target)
guide_targets <- guides %>% 
  mutate(revmap = seq_len(nrow(guides))) %>% 
  left_join(targets_tbl, by = "revmap") %>% 
  select(-c(score, target_score, revmap))

# write guide targets to output file
write_tsv(guide_targets, file = snakemake@output$guide_targets)

# create bed file containing target annotations
write_tsv(distinct(select(targets_tbl, -revmap)), file = snakemake@output$targets_bed,
          col_names = FALSE)

## Merge guide targets from both TAP-seq experiments and write to new files

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# column types in guide targets files
guide_targets_cols <- cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  strand = col_character(),
  spacer = col_character(),
  target_chr = col_character(),
  target_start = col_integer(),
  target_end = col_integer(),
  target_name = col_character(),
  target_strand = col_character(),
  target_type = col_character()
)

# load guide target files for both TAP-seq experiments
guide_targets_files <- c(chr8 = snakemake@input$chr8, chr11 = snakemake@input$chr11)
guide_targets <- lapply(guide_targets_files, FUN = read_tsv, col_types = guide_targets_cols,
                        progress = FALSE)

# combine into one table and remove any duplicates due to controls used in both experiments
guide_targets <- distinct(bind_rows(guide_targets))

# create bed format guide coordinates
guides_bed <- guide_targets %>% 
  mutate(score = 500) %>% 
  select(chr, start, end, name, score, strand, spacer) %>% 
  arrange(chr, start, end)

# create bed format target coordinates
targets_bed <- guide_targets %>% 
  mutate(score = 500) %>% 
  select(target_chr, target_start, target_end, target_name, score, target_strand) %>% 
  distinct() %>% 
  arrange(target_chr, target_start, target_end)

# write output to files
write_tsv(guide_targets, file = snakemake@output$guide_targets)
write_tsv(guides_bed, file = snakemake@output$guides_bed, col_names = FALSE)
write_tsv(targets_bed, file = snakemake@output$targets_bed, col_names = FALSE)

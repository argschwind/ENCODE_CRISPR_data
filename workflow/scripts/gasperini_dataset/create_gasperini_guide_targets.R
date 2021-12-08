## Create file with target annotations for each gRNA based on targets defined by Gasperini et al.

# required packages
library(dplyr)
library(readxl)
library(readr)

# column types in guides file
guides_cols <- cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character(),
  spacer = col_character()
)

# load re-mapped guides
guides <- read_tsv(snakemake@input$guides, col_types = guides_cols,
                   col_names = names(guides_cols$cols))

# load table containing at scale screen guides and their targets
targets <- read_excel(snakemake@input$targets, sheet = "S2A_AtScale_library_gRNA.cs")

# extract relevant columns from targets containing guide target coordinates
targets_coords <- targets %>%
  select(spacer = Spacer, target_chr = chr.candidate_enhancer,
         target_start = start.candidate_enhancer, target_end = stop.candidate_enhancer,
         target_name = Target_Site) %>% 
  mutate(target_strand = "*")

# add targets to guides based on spacer sequence to create output targets file
guide_targets <- guides %>% 
  select(-strand) %>% 
  left_join(targets_coords, by = "spacer")

# create arbitrary target regions for positive controls and TSS (guides that were not designed
# against putative enhancers)
guide_targets <- guide_targets %>%
  group_by(target_name) %>%
  mutate(target_chr = if_else(is.na(target_chr), true = chr, false = target_chr),
         target_start = if_else(is.na(target_start), true = min(start),
                                false = as.integer(target_start)),
         target_end = if_else(is.na(target_end), true = max(end), false = as.integer(target_end)))

# extract Gasperini target coordinates in bed format
gasperini_targets <- guide_targets %>% 
  mutate(target_score = 666) %>% 
  select(target_chr, target_start, target_end, target_name, target_score, target_strand) %>% 
  distinct() %>% 
  arrange(target_chr, target_start)

# save guide targets to file
write_tsv(guide_targets, file = snakemake@output$guide_targets)

# write Gasperini targets to bed file
write_tsv(gasperini_targets, file = snakemake@output$targets_bed, col_names = FALSE)

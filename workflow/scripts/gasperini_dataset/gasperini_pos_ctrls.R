## Extract Gasperini positive TSS control targets and guides, and save to file

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

# Load input data ----------------------------------------------------------------------------------

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

# load supplementary data table with Gasperini gRNA library including targets
targets <- read_excel(snakemake@input$targets, sheet = "S2A_AtScale_library_gRNA.cs")

# add G to 19 bp long spacer sequences in targets table
targets <- targets %>% 
  mutate(Spacer = if_else(nchar(Spacer) == 19, true = paste0(Spacer, "G"), false = Spacer))

# Extract positive control targets and guides ------------------------------------------------------

# get positive TSS targeting controls
pos_ctrls <- targets %>% 
  filter(Category %in% "TSS" | grepl(Target_Site, pattern = "tss")) %>% 
  select(Spacer, Target_Site, Category)

# add re-mapped guide RNA positions
pos_ctrls <- inner_join(pos_ctrls, guides, by = c("Spacer" = "spacer"))

# create coordinates for every control target based on guide coordinates
pos_ctrl_targets <- pos_ctrls %>% 
  dplyr::rename(target_chr = chr) %>% 
  group_by(Target_Site, Category, target_chr) %>% 
  summarize(target_start = min(start), target_end = max(end), .groups = "drop")

# add target coordinates to positive controls
pos_ctrls <- left_join(pos_ctrls, pos_ctrl_targets, by = c("Target_Site", "Category"))

# reformat to the same format as guide targets file
pos_ctrls <- pos_ctrls %>% 
  dplyr::rename(spacer = Spacer) %>% 
  mutate(target_name = paste("TSSCtrl", Target_Site, sep = "|"),
         target_strand = ".", target_type = "TSSCtrl") %>% 
  select(name, target_chr, target_start, target_end, target_name, target_strand, target_type)

# write new guide targets to file
write_tsv(pos_ctrls, file = snakemake@output[[1]])

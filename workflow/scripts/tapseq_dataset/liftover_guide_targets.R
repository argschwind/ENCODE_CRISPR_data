## Lift guide targets file from hg19 to hg38 using lifted over guide and target coordinates

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# column types in guide targets file
guide_targets_cols <- cols(
  .default = col_character(),  # for target_type columns, which is optional as of now
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
  target_strand = col_character()
)

# load hg19 guide targets file
guide_targets_hg19 <- read_tsv(snakemake@input$guide_targets_hg19, col_types = guide_targets_cols)

# original column names need later to format output
guide_targets_colnames <- colnames(guide_targets_hg19)

# column to load first 6 columns of a bed format file (ignores any additional columns)
bed_cols <- cols_only(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character()
)

# load lifted over guide coordinates
guides_hg38 <- read_tsv(snakemake@input$guides_hg38, col_types = bed_cols,
                        col_names = names(bed_cols$cols))

# load lifted over target coordinates
targets_hg38 <- read_tsv(snakemake@input$targets_hg38, col_types = bed_cols,
                         col_names = names(bed_cols$cols))

# rename columns in targets to match guide targets file
colnames(targets_hg38) <- paste0("target_", colnames(targets_hg38))

# replace hg19 guide coordinates with hg38 coordinates
guide_targets_hg38 <- guide_targets_hg19 %>% 
  select(-c(chr, start, end, strand)) %>% 
  left_join(select(guides_hg38, -score), by = "name")

# replace hg19 target coordinates with hg38 coordinates
guide_targets_hg38 <- guide_targets_hg38 %>% 
  select(-c(target_chr, target_start, target_end, target_strand)) %>% 
  left_join(select(targets_hg38, -target_score), by = "target_name")

# remove any NAs from coordinates that could not be lifted over and rearrange columns
output <- guide_targets_hg38 %>% 
  drop_na(chr, start, end, target_chr, target_start, target_end) %>% 
  select(all_of(colnames(guide_targets_hg19)))

# raise warning if some guides or targets were not lifted over
if (nrow(output) < nrow(guide_targets_hg38)) {
  warning(nrow(guide_targets_hg38) - nrow(output),
          "entries removed due to guides or targets not being lifted over", call. = FALSE)
}

# write output to file
write_tsv(output, file = snakemake@output[[1]])

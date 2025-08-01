## Process ENCODE portal metadata file to extract download URLs

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# load metadata file
meta <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# filter for element-level quantification files and extract relevant columns
output <- meta %>% 
  filter(`File format` == "bed CRISPR element quantifications") %>% 
  select(accession = `File accession`, url = `File download URL`)

# save to output file
write_tsv(output, snakemake@output[[1]])

## Filter EP benchmarking dataset based on minimum statistical power for specified effect size. Adds
## a low power flag to the ValidConnection column if a pair is a negative below a specified power
## threshold (positives pass filter regardless of power). If specified removes all invalid
## connections based on ValidConnection column before writing to output.

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# column types in EP benchmarking datasets
crispr_data_cols <- cols(
  .default = col_double(),  # for power columns, which can vary in names
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
  ValidConnection = col_character(),
  CellType = col_character(),
  Reference = col_character(),
  Regulated = col_logical()
)

# load unfiltered dataset
dat <- read_tsv(snakemake@input[[1]], col_types = crispr_data_cols, progress = FALSE)

# power filter
power_filt_col <- paste0("PowerAtEffectSize", snakemake@wildcards$es)
min_power <- as.numeric(snakemake@wildcards$pwr)

# ValidConnection flag to add to low power pairs
power_flag <- paste0("< ", min_power * 100, "% power at ", snakemake@wildcards$es, "% effect size")

# add filter data based on power
output <- dat %>% 
  mutate(ValidConnection = case_when(
    ValidConnection != "TRUE" ~ ValidConnection,
    (!!sym(power_filt_col) < 0.8 | is.na(!!sym(power_filt_col))) & Significant != TRUE ~ power_flag,
    TRUE ~ "TRUE"
  ))

# filter out any invalid connections if specified
if (snakemake@params$remove_filtered_pairs == TRUE) {
  output <- filter(output, ValidConnection == "TRUE")
}

# write output to file
write_tsv(output, file = snakemake@output[[1]])

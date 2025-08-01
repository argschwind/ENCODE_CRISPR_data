## Filter EP benchmarking datasets based on expected indirect effects beyond a given distance

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# distance threshold for filtering
dist_threshold <- as.numeric(snakemake@wildcards$dist) * 1000

crispr <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# extract original column names for later formatting
crispr_cols <- colnames(crispr)

# calculate distance to TSS
crispr <- crispr %>%
  mutate(distToTSS = if_else(
    chrom == chrTSS,
    true = abs(((chromStart + chromEnd) / 2) - ((startTSS + endTSS) / 2)),
    false = NA_real_))

# Hard filter for any E-G pairs over distance threshold --------------------------------------------

# filter crispr data based on distance threshold
crispr_filt <- crispr %>% 
  filter(distToTSS <= dist_threshold) %>% 
  select(all_of(crispr_cols))

# save to output file in EPbenchmarking format for CRISPR benchmarking
write_tsv(crispr_filt, file = snakemake@output$filt)

# Re-assign positives over distance threshold as negatives -----------------------------------------

# set 'Regulated' to FALSE for all pairs over distance threshold 
crispr_flip <- crispr %>% 
  mutate(Regulated = if_else(distToTSS > dist_threshold, true = FALSE, false = Regulated)) %>% 
  select(all_of(crispr_cols))

# save to output file
write_tsv(crispr_flip, file = snakemake@output$flip)

# Filter out positives plus a proportionate fraction of negatives over distance threshold ----------

# get all positive and negatives
positives <- filter(crispr, Regulated == TRUE)
negatives <- filter(crispr, Regulated == FALSE)

# filter positives based on distance threshold and get proportion of filtered pairs
positives_filt <- filter(positives, distToTSS <= dist_threshold)
filt_prop <- nrow(positives_filt) / nrow(positives)

# randomly sample an equal proportion of negatives
negatives_sample <- sample_n(negatives, size = round(nrow(negatives) * filt_prop))

# combine filtered positives and sampled negatives to form output
crispr_prop_filt <- bind_rows(positives_filt, negatives_sample) %>% 
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol) %>% 
  select(all_of(crispr_cols))

# save to output file
write_tsv(crispr_prop_filt, file = snakemake@output$prop_filt)

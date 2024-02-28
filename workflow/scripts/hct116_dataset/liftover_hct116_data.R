
# save.image("liftover_hct116.rda")
# stop()

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

# load results
results <- read_tsv(snakemake@input$results, show_col_types = FALSE)

# load lifted over enhancers
enh_cols <- c("chrom", "chromStart", "chromEnd", "PerturbationTargetID", "score", "strand")
enh <- read_tsv(snakemake@input$enh, col_names = enh_cols, show_col_types = FALSE)

# replace hg19 enhancer coordinates with hg38 coordinates
results <- results %>% 
  select(-c(chrom, chromStart, chromEnd)) %>% 
  left_join(select(enh, -c(score, strand)), by = "PerturbationTargetID") %>% 
  relocate(chrom, chromStart, chromEnd, .before = "name") %>% 
  filter(!is.na(chrom))

# load TSS annotations
tss_cols <- c("chrTSS", "startTSS", "endTSS", "measuredGeneSymbol", "score", "strandGene")
tss <- read_tsv(snakemake@input$tss, col_names = tss_cols, show_col_types = FALSE)

# replace hg19 enhancer coordinates with hg38 coordinates
results <- results %>% 
  left_join(select(tss, "chrTSS", "startTSS", "endTSS", "strandGene", "measuredGeneSymbol"),
            by = "measuredGeneSymbol") %>% 
  relocate(chrTSS, startTSS, endTSS, strandGene, .after = "PerturbationTargetID")

# compute distance to TSS
results <- results %>% 
  mutate(enh_center = (chromStart + chromEnd) / 2) %>% 
  mutate(distToTSS = abs(enh_center - startTSS))

# Create EPbenchmarking output ---------------------------------------------------------------------

# filter based on distance to TSS
dist_filter <- c(1000, 1e6)
results <- filter(results, abs(distToTSS) > dist_filter[[1]], abs(distToTSS) <= dist_filter[[2]])

# select all columns for EPbenchmarking output
results <- results %>% 
  mutate(CellType = "HCT116") %>% 
  select(chrom, chromStart, chromEnd, name, EffectSize, chrTSS, startTSS,
         endTSS, measuredGeneSymbol, Significant, pValueAdjusted, starts_with("PowerAtEffectSize"),
         ValidConnection, CellType, Reference, Regulated) %>% 
    arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)

write_tsv(results, file = snakemake@output[[1]])

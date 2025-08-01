## Calculate cis and trans positive hit ratios for a sample. Cis positive ratio is calculated across
## specified distance bins

# save.image("calc_rate.rda")
# stop()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# load cis and trans analysis results
trans_results <- read_tsv(snakemake@input$trans_results, show_col_types = FALSE)
cis_results <- read_tsv(snakemake@input$cis_results, show_col_types = FALSE)

# load processed ENCODE format results including filters
encode_results <- read_tsv(snakemake@input$encode_results, show_col_types = FALSE)

# get gene symbol for every gene id from processed ENCODE data
genes <- encode_results %>%
  select(gene = measuredEnsemblID, gene_symbol = measuredGeneSymbol) %>% 
  distinct()

# add gene symbols to cis results 
cis_results <- left_join(cis_results, genes, by = "gene")

# add ValidConnection column to cis results
cis_results <- cis_results %>% 
  mutate(name = paste0(gene_symbol, "|", pert_chr, ":", pert_start, "-", pert_end, ":.")) %>% 
  left_join(distinct(select(encode_results, name, ValidConnection)), by = "name")

# get valid enhancers only
enh_filter <- c("overlaps potential promoter", "TSS targeting guide(s)")
valid_enh <- encode_results %>% 
  filter(!ValidConnection %in% enh_filter) %>% 
  pull(PerturbationTargetID)

# add ValidConnection column to trans-results
trans_results <- trans_results %>% 
  mutate(pert_id = paste0(pert_chr, ":", pert_start, "-", pert_end, ":.")) %>% 
  mutate(ValidConnection = if_else(pert_id %in% valid_enh, true = "TRUE", false = "Invalid enhancer"))

# filter cis and trans results based on ValidConnection
cis_results_filt <- filter(cis_results, ValidConnection == "TRUE")
trans_results_filt <- filter(trans_results, ValidConnection == "TRUE")

# add significance and regulated columns to cis results
cis_results_filt <- cis_results_filt %>% 
  mutate(significant = pval_adj < 0.05,
         regulated_negative = significant & logFC < 0,
         regulated_positive = significant & logFC >= 0) %>% 
  filter(!is.na(regulated_negative), !is.na(regulated_positive))

# distance bins for calculating positive hit ratio
cis_results_filt <- mutate(cis_results_filt, abs_dist_to_tss = abs(dist_to_tss))
max_dist <- ceiling(max(cis_results_filt$abs_dist_to_tss, na.rm = TRUE) / 1e6) * 1e6
dist_bins <- seq(0, max_dist, by = snakemake@params$bin_size)

# bin cis pairs by distance
cis_results_filt <- cis_results_filt %>% 
  mutate(dist_bin = cut(abs_dist_to_tss, breaks = dist_bins, include.lowest = TRUE))

# calculate positive ratios and mean distance per distance bin
cis_pos_ratio <- cis_results_filt %>% 
  group_by(dist_bin) %>% 
  summarize(total_pairs = n(),
            significant_pairs = sum(significant, na.rm = TRUE),
            negative_pairs = sum(regulated_negative, na.rm = TRUE),
            positive_pairs = sum(regulated_positive, na.rm = TRUE),
            positive_rate_significant = significant_pairs / total_pairs,
            positive_rate_negative = negative_pairs / total_pairs,
            positive_rate_positive = positive_pairs / total_pairs,
            mean_dist = mean(abs_dist_to_tss))

# calculate trans-effects positive ratio
trans_pos_ratio <- trans_results_filt %>% 
  summarize(total_pairs = n(),
            significant_pairs = sum(significant, na.rm = TRUE),
            negative_pairs = sum(regulated_negative, na.rm = TRUE),
            positive_pairs = sum(regulated_positive, na.rm = TRUE),
            positive_rate_significant = significant_pairs / total_pairs,
            positive_rate_negative = negative_pairs / total_pairs,
            positive_rate_positive = positive_pairs / total_pairs)

# save output to files
write_tsv(cis_pos_ratio, snakemake@output$cis_ratio)
write_tsv(trans_pos_ratio, snakemake@output$trans_ratio)

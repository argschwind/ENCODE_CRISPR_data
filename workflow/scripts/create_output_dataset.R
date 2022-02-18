## Combine differential expression and power simulation results to create output data sets

# save.image("create_dataset.rda")
# stop()

# required packages
library(tidyverse)
library(rtracklayer)

# load input data ----------------------------------------------------------------------------------

# column types in DE results 
diff_expr_cols <-cols(
  perturbation = col_character(),
  gene = col_character(),
  logFC = col_double(),
  ci_high = col_double(),
  ci_low = col_double(),
  pvalue = col_double(),
  pval_adj = col_double(),
  cells = col_double(),
  avg_expr = col_double(),
  chr = col_character(),
  pert_start = col_double(),
  pert_end = col_double(),
  gene_tss = col_double(),
  gene_strand = col_character(),
  dist_to_tss = col_double()
)

# load differential expression results
diff_expr <- read_csv(snakemake@input$diff_expr, col_types = diff_expr_cols, progress = FALSE)

# column types in power files
power_cols <- cols(
  perturbation = col_character(),
  gene = col_character(),
  disp_outlier_deseq2 = col_logical(),
  power = col_double(),
  effect_size = col_double()
)

# load power files and combine into one data frame
power <- snakemake@input$power_sim %>% 
  lapply(FUN = read_csv, col_types = power_cols, progress = FALSE) %>% 
  bind_rows()

# load genome annotations
annot <- import(snakemake@input$annot)

# remove any version numbers from gene ids
annot$gene_id <- sub("\\..+$", "", annot$gene_id)

## HACK: remove wrong hemoglobon annotations, since they are too long
hemoglobins <- c("HBE1", "HBG2")
annot <- annot[!annot$gene_name %in% hemoglobins]

# column types in guide targets file
guide_targets_cols <- cols(
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  name = col_character(),
  strand = col_character(),
  spacer = col_character(),
  target_chr = col_character(),
  target_start = col_double(),
  target_end = col_double(),
  target_name = col_character(),
  target_strand = col_character()
)

# load guide targets file
guide_targets <- read_tsv(snakemake@input$guide_targets, col_types = guide_targets_cols,
                          progress = FALSE)

# add simulated power to differential expression results -------------------------------------------

# reformat effect size to percentage and transform to wide format with new column names
power <- power %>%
  mutate(effect_size = effect_size * 100) %>% 
  pivot_wider(names_from = effect_size, names_prefix = "PowerAtEffectSize", values_from = power)

# add power to DE results
output <- left_join(diff_expr, select(power, -disp_outlier_deseq2), by = c("perturbation", "gene"))

# add gene names respectively ensembl ids to results -----------------------------------------------

# create gene name - ensembl id table
gene_ids <- mcols(annot)[, c("gene_id", "gene_name")] %>% 
  as.data.frame() %>% 
  distinct() %>% 
  dplyr::rename(measuredEnsemblID = gene_id, measuredGeneSymbol = gene_name)

# add missing gene id to data and rename columns according to ENCODE format
if (snakemake@params$gene_ids == "ensembl") {
  
  output <- output %>% 
    dplyr::rename(measuredEnsemblID = gene) %>% 
    left_join(gene_ids, "measuredEnsemblID")
  
} else if (snakemake@params$gene_ids == "symbol") {
  
  output <- output %>% 
    dplyr::rename(measuredGeneSymbol = gene) %>% 
    left_join(gene_ids, "measuredGeneSymbol")
  
} else {
  stop("Invalid 'gene_ids' snakemake parameter.", call. = FALSE)
}

# identify and tag promoter targeting guides, and intronic and promoter enhancers ------------------

# create GRangesList with all annotated exons per annotated gene and transcripts
exons <- annot[annot$type == "exon"]
genes <- split(exons, f = exons$gene_id)
txs <- split(exons, f = exons$transcript_id)

# get gene locus coordinates
gene_loci <- unlist(range(genes))

# get all promoters for annotated transcripts
promoters <- unlist(promoters(range(txs), upstream = snakemake@params$tss_min_dist))

# get gene id for every transcript id
gene_ids_per_tx_id <- mcols(exons)[, c("transcript_id", "gene_id")] %>% 
  as.data.frame() %>% 
  distinct() %>% 
  deframe()

# add gene id to promoters metadata columns
mcols(promoters) <- DataFrame(feature_id = gene_ids_per_tx_id[names(promoters)], type = "promoter")

# add metadata columns to exons and gene loci
mcols(exons) <- DataFrame(feature_id = exons$gene_id, type = "exon")
mcols(gene_loci) <- DataFrame(feature_id = names(gene_loci), type = "gene")

# combine exons, genes and promoters into one feature GRanges object
features <- c(exons, gene_loci, promoters)
names(features) <- NULL

# extract perturbation coordinates
pert_coords <- output %>% 
  select(chr, pert_start, pert_end, perturbation) %>% 
  distinct() %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

# overlap perturbations with features
feature_ovl <- findOverlapPairs(pert_coords, features, ignore.strand = TRUE)

# create table with overlapping features per perturbation in long format
feature_ovl <- tibble(perturbation = first(feature_ovl)$perturbation,
                      feature_id = second(feature_ovl)$feature_id,
                      type = second(feature_ovl)$type)

# add tested genes for perturbations overlapping any features
feature_ovl <- feature_ovl %>% 
  distinct() %>% 
  left_join(select(output, perturbation, measuredEnsemblID), by = "perturbation")

# summarize feature overlaps
feature_ovl_summary <- feature_ovl %>% 
  group_by(perturbation, measuredEnsemblID, type) %>% 
  summarize(ovl = n() > 0,
            ovl_target = unique(measuredEnsemblID) %in% feature_id,
            .groups = "drop")

# reformat overlap summary
feature_ovl_summary <- feature_ovl_summary %>% 
  pivot_longer(cols = c(ovl, ovl_target), names_to = "ovl_type", values_to = "ovl") %>% 
  unite(type, ovl_type, type, sep = "_") %>% 
  pivot_wider(names_from = type, values_from = ovl, values_fill = FALSE)

# infer whether perturbation overlaps intron
feature_ovl_summary <- feature_ovl_summary %>% 
  mutate(ovl_intron = if_else((ovl_gene & !ovl_exon), true = TRUE, false = FALSE),
         ovl_target_intron = if_else((ovl_intron & ovl_target_gene), true = TRUE, false = FALSE))

# add selected overlaps to output
output <- feature_ovl_summary %>% 
  select(perturbation, measuredEnsemblID, ovl_target_exon, ovl_target_intron) %>% 
  left_join(output, ., by = c("perturbation", "measuredEnsemblID"))

# create ValidConnection column based on overlaps
output <- output %>% 
  mutate(ValidConnection = case_when(
    ovl_target_exon == TRUE ~ "overlaps target gene exon",
    ovl_target_intron == TRUE ~ "overlaps target gene intron",
    TRUE ~ "TRUE"
  ))

# add TSS targeting flag to known TSS controls (if TSS control pattern is provided)
if (!is.null(snakemake@params$tss_ctrl_tag)) {
  output <- output %>% 
    mutate(ValidConnection = if_else(grepl(perturbation, pattern = snakemake@params$tss_ctrl_tag),
                                     true = "TSS targeting guide(s)", false = ValidConnection))
}

# add gRNA sequences -------------------------------------------------------------------------------

# get whether perturbation id refers to guide or target name
if (snakemake@wildcards$strategy == "perGRNA") {
  pert_id <- "name"
} else {
  pert_id <- "target_name"
}

# get all gRNA sequences per perturbation and concatenate to one string
grna_seqs <- guide_targets %>% 
  select(perturbation = all_of(pert_id), spacer) %>% 
  group_by(perturbation) %>% 
  summarize(guideSeq = paste(spacer, collapse = ";"))

# add to data
output <- left_join(output, grna_seqs, by = "perturbation")

# add NA for 'guideSpacerSeq'
output <- mutate(output, guideSpacerSeq = NA_character_)

# create ENCODE style output -----------------------------------------------------------------------

# add significance column and gene tss coordinates
output <- output %>% 
  mutate(Significant = pval_adj < snakemake@params$padj_threshold, 
         startTSS = gene_tss, endTSS = gene_tss + 1, strandGene = gene_strand)

# add additional meta data columns in ENCODE style to create output
output <- output %>% 
  mutate(strandPerturbationTarget = ".",
         PerturbationTargetID = paste0(chr, ":", pert_start, "-", pert_end, ":",
                                       strandPerturbationTarget),
         name = paste(measuredGeneSymbol, PerturbationTargetID, sep = "|"),
         CellType = snakemake@params$cell_type,
         Reference = snakemake@params$reference)

# get columns containing simulated power
power_colnames <- grep(colnames(output), pattern = "^PowerAtEffectSize.+", value = TRUE)

# reformat to ENCODE style (last 3 columns are not ENCODE format and need to be stripped for upload)
output <- output %>% 
  select(chrom = chr, chromStart = pert_start, chromEnd = pert_end, name, EffectSize = logFC,
         strandPerturbationTarget, PerturbationTargetID, chrTSS = chr, startTSS,
         endTSS, strandGene, EffectSize95ConfidenceIntervalLow = ci_low,
         EffectSize95ConfidenceIntervalHigh = ci_high, measuredGeneSymbol, measuredEnsemblID,
         guideSpacerSeq, guideSeq, Significant, pValue = pvalue, pValueAdjusted = pval_adj,
         all_of(power_colnames), ValidConnection, CellType, Reference, distToTSS = dist_to_tss,
         avgGeneExpr = avg_expr, nPertCells = cells)

# make sure output is sorted according to genomic coordinates
output <- arrange(output, chrom, chromStart, chromEnd, measuredGeneSymbol)

# write to .tsv file
write_tsv(output, file = snakemake@output[[1]])

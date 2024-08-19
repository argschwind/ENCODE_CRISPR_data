## Reformat processed HCT116 FlowFISH results to ENCODE CRISPR format

# required packages
suppressWarnings({
  library(tidyverse)
  library(rtracklayer)
})


# Process results and add all required columns -----------------------------------------------------

# load processed HCT116 results
positives <- read_csv(snakemake@input$pos, show_col_types = FALSE)
negatives <- read_csv(snakemake@input$neg, show_col_types = FALSE)

# combine into one table and add Regulated column
results <- bind_rows(pos = positives, neg = negatives, .id = "Regulated") %>% 
  mutate(Regulated = if_else(Regulated == "pos", true = TRUE, false = FALSE))

# load TSS annotations
tss_cols <- c("chrTSS", "startTSS", "endTSS", "TargetGene", "score", "strandGene")
tss <- read_tsv(snakemake@input$tss, col_names = tss_cols, show_col_types = FALSE)

# load genome annotations and remove any version numbers from gene ids
annot <- import(snakemake@input$annot)
annot$gene_id <- sub("\\..+$", "", annot$gene_id)

# create gene name - ensembl id table
gene_ids <- mcols(annot)[, c("gene_id", "gene_name")] %>% 
  as.data.frame() %>% 
  distinct()

# manually change the name for 'ITPRID2', which is called 'SSFA2' in results
gene_ids <- mutate(gene_ids, gene_name = if_else(gene_name == "ITPRID2", "SSFA2", gene_name))

# add gene ids to results
results <- left_join(results, gene_ids, by = c("TargetGene" = "gene_name"))

# add TSS coordinates to results
results <- results %>% 
  left_join(select(tss, "chrTSS", "startTSS", "endTSS", "strandGene", "TargetGene"),
            by = "TargetGene")

# calculate percent change effect size
results <- results %>% 
  mutate(EffectSize = -EnhancerEffect_noAux)

# label significant E-G pairs
sig_pval <- snakemake@params$padj_threshold
results <- results %>% 
  mutate(Significant = (FDR_pval_noAux_rep1 < sig_pval & FDR_pval_noAux_rep2 < sig_pval))


# Reformat into EP benchmarking format -------------------------------------------------------------

# Create other required columns with default values
results <- results %>% 
  mutate(strandPerturbationTarget = ".",
         name = paste0(TargetGene, "|", name_hg19),
         EffectSize95ConfidenceIntervalLow = NA_real_,
         EffectSize95ConfidenceIntervalHigh = NA_real_,
         measuredEnsemblID = NA_character_,
         guideSpacerSeq = NA_character_,
         guideSeq = NA_character_,
         ValidConnection = "TRUE",
         Notes = NA_character_,
         Reference = snakemake@params$reference)

# select ENCODE format columns
output <- results %>% 
  select(chrom = chr_hg38, chromStart = start_hg38, chromEnd = end_hg38, name, EffectSize,
         strandPerturbationTarget, PerturbationTargetID = name_hg19, chrTSS, startTSS, endTSS,
         strandGene, EffectSize95ConfidenceIntervalLow, EffectSize95ConfidenceIntervalHigh,
         measuredGeneSymbol = TargetGene, measuredEnsemblID, guideSpacerSeq, guideSeq, Significant,
         pValue = pval_noAux, pValueAdjusted = FDR_pval_noAux,
         PowerAtEffectSize25 = `Power_25%_effect`, ValidConnection, Notes, Reference,
         distToTSS = distance, Regulated)

# make sure output is sorted according to genomic coordinates
output <- arrange(output, chrom, chromStart, chromEnd, measuredGeneSymbol)

# write to .tsv file
write_tsv(output, file = snakemake@output[[1]])

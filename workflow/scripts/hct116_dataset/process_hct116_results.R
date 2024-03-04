## Reformat processed HCT116 FlowFISH results to ENCODE CRISPR format


# required packages
suppressWarnings({
  library(tidyverse)
  library(rtracklayer)
})


# Process results and add all required columns -----------------------------------------------------

# load processed HCT116 results
results <- read_csv(snakemake@input$results, show_col_types = FALSE)

# filter out any negative control perturbations
results <- filter(results, Category != "negative_control")

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

# compute distance to TSS
results <- results %>% 
  mutate(enh_center_hg38 = (start_hg38 + end_hg38) / 2) %>% 
  mutate(distToTSS = round(abs(enh_center_hg38 - startTSS)))

# calculate percent change effect size
results <- results %>% 
  mutate(EffectSize = -EnhancerEffect_noAux)

# label significant E-G pairs
results <- mutate(results, Significant = FDR_pval_noAux < snakemake@params$padj_threshold)


# Identify enhancers in target genes and potential promoters ---------------------------------------

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
pert_coords <- results %>% 
  select(chr = chr_hg38, start = start_hg38, end = end_hg38, perturbation = name_hg19) %>% 
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
  left_join(select(results, name_hg19, gene_id), by = c("perturbation" = "name_hg19"))

# summarize feature overlaps
feature_ovl_summary <- feature_ovl %>% 
  group_by(perturbation, gene_id, type) %>% 
  summarize(ovl = n() > 0,
            ovl_target = unique(gene_id) %in% feature_id,
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

# add selected overlaps to results
results <- feature_ovl_summary %>% 
  select(perturbation, gene_id, ovl_target_exon, ovl_target_intron, ovl_promoter) %>% 
  left_join(results, ., by = c("name_hg19" = "perturbation", "gene_id"))

# create ValidConnection column based on overlaps
results <- results %>% 
  mutate(ValidConnection = case_when(
    ovl_target_exon == TRUE ~ "overlaps target gene exon",
    ovl_target_intron == TRUE ~ "overlaps target gene intron",
    ovl_promoter == TRUE ~ "overlaps potential promoter",
    TRUE ~ "TRUE"
  ))


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
         Notes = NA_character_,
         Reference = snakemake@params$reference)

# select ENCODE format columns
output <- results %>% 
  select(chrom = chr_hg38, chromStart = start_hg38, chromEnd = end_hg38, name, EffectSize,
         strandPerturbationTarget, PerturbationTargetID = name_hg19, chrTSS, startTSS, endTSS,
         strandGene, EffectSize95ConfidenceIntervalLow, EffectSize95ConfidenceIntervalHigh,
         measuredGeneSymbol = TargetGene, measuredEnsemblID, guideSpacerSeq, guideSeq, Significant,
         pValue = pval_noAux, pValueAdjusted = FDR_pval_noAux,
         PowerAtEffectSize25 = `Power_25%_effect`, ValidConnection, Notes, Reference, distToTSS)

# make sure output is sorted according to genomic coordinates
output <- arrange(output, chrom, chromStart, chromEnd, measuredGeneSymbol)

# write to .tsv file
write_tsv(output, file = snakemake@output[[1]])

## Download HCR-FlowFish files from ENCODE portal and combine into one file

# save.image("dl.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(GenomicRanges)
  library(rtracklayer)
  library(BiocParallel)
})

## Download HCR-FlowFISH data ----------------------------------------------------------------------

# register parallel backend if specified (if more than 1 thread provided)
if (snakemake@threads > 1) {
  message("Registering parallel backend with ", snakemake@threads, " cores.")
  register(MulticoreParam(workers = snakemake@threads))
} else {
  message("Registering serial backend.")
  register(SerialParam())
}

# load file describing ENCODE CRISPR format
format <- read_tsv(snakemake@input$encode_format, show_col_types = FALSE)

# load parsed file with download urls
urls <- read_tsv(snakemake@input$urls,
                 col_types = cols(accession = col_character(), url = col_character()))

# create vector with download urls
urls <- deframe(urls)

# download and read all files and combine into on data frame
dl <- urls %>% 
  bplapply(FUN = read_tsv, col_names = FALSE, show_col_types = FALSE) %>% 
  bind_rows()

# get minimal column names for ENCODE element-level format
min_cols <- format %>% 
  filter(`Required for cCRE-gene analysis (CRISPR element_quantifications BED3+22 format)` == TRUE) %>%
  pull(ColumnName)

# set these as column names
colnames(dl) <- min_cols


## Reformat to align with other CRISPR datasets ----------------------------------------------------

# make sure that NAs are NA and not "."
dl <- dl %>% 
  mutate(guideSpacerSeq = if_else(guideSpacerSeq == ".", true = NA_character_, false = guideSpacerSeq),
         guideSeq = if_else(guideSeq == ".", true = NA_character_, false = guideSeq)) %>% 
  mutate(EffectSize = abs(EffectSize) * -1)

# add gene symbol to 'name'
dl <- unite(dl, col = name, measuredGeneSymbol, name, sep = "|", remove = FALSE)

# remove version number from gene ids
dl <- mutate(dl, measuredEnsemblID = sub("\\..+", "", measuredEnsemblID))

# load genome annotations, filter out ignored transcripts and remove version numbers from gene ids
annot <- import(snakemake@input$annot)
annot <- annot[!annot$transcript_id %in% snakemake@params$ignore_txs]
annot$gene_id <- sub("\\..+$", "", annot$gene_id)

# get TSS annotations for all target genes
exons <- annot[annot$type == "exon"]
genes <- split(exons, f = exons$gene_name)
tss <- promoters(unlist(range(genes)), upstream = 0, downstream = 1)
  
# convert to data frame
tss <- tibble(chrTSS = as.character(seqnames(tss)), startTSS = start(tss) - 1, endTSS = end(tss),
              measuredGeneSymbol = names(tss)) %>% 
  distinct()

# add to HCR-FlowFISH data
dl <- dl %>% 
  select(-c(chrTSS, startTSS, endTSS)) %>% 
  left_join(tss, by = "measuredGeneSymbol")

# add reference and distance to TSS (latter not ENCODE standard format) columns
dl <- dl %>% 
  mutate(ValidConnection = as.character(ValidConnection)) %>% 
  mutate(Reference = snakemake@params$reference) %>% 
  mutate(pert_center = floor((chromStart + chromEnd) / 2),
         distToTSS = case_when(
           chrom != chrTSS ~ NA_real_,
           strandGene == "-" ~ case_when(
             pert_center < startTSS ~ startTSS - pert_center,
             pert_center > endTSS ~ endTSS - pert_center,
             TRUE ~ 0),
           strandGene == "+" ~ case_when(
             pert_center < startTSS ~ pert_center - startTSS,
             pert_center > endTSS ~ pert_center - endTSS,
             TRUE ~ 0),
           TRUE ~ NA_real_
         )) %>% 
  select(-pert_center)


## Filter for E-G pairs not overlapping E-G pairs in the original combined dataset -----------------


# load combined CRISPR dataset
combined <- read_tsv(snakemake@input$training_crispr_data, show_col_types = FALSE)

# create GRanges objects for E-G pairs in the HCR-FlowFISH and combined dataset
hcr_gr <- with(dl,
               GRanges(seqnames = paste0(chrom,":", measuredGeneSymbol),
                       ranges = IRanges(chromStart, chromEnd), name = name))
combined_gr <- with(combined,
                    GRanges(seqnames = paste0(chrom,":", measuredGeneSymbol),
                            ranges = IRanges(chromStart, chromEnd)))

# set common seqlevels to prevent warnings during overlap
seqlevels_all_pairs <- as.character(unique(c(seqnames(hcr_gr), seqnames(combined_gr))))
seqlevels(hcr_gr) <- seqlevels_all_pairs
seqlevels(combined_gr) <- seqlevels_all_pairs

# find HCR-FlowFISH E-G pairs that overlap combined CRISPR dataset
ovl <- subsetByOverlaps(hcr_gr, combined_gr)

# add filter to these in ValidConnection column
dl <- dl %>% 
  mutate(ValidConnection = if_else(name %in% ovl$name, true = "Overlaps training data",
                                   false = ValidConnection))


## Filter for elements overlapping candidate elements only -----------------------------------------

# load candidate elements
cres <- import(snakemake@input$cres)

# create GRanges object for all perturbed HCR-FlowFISH elements
pert_coords <- dl %>% 
  select(chrom, chromStart, chromEnd, PerturbationTargetID) %>%
  distinct() %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

# find HCR-FlowFISH elements that overlap combined CRISPR dataset
ovl <- subsetByOverlaps(pert_coords, cres)

# add filter to these in ValidConnection column
dl <- dl %>% 
  mutate(ValidConnection = if_else(PerturbationTargetID %in% ovl$PerturbationTargetID,
                                   true = ValidConnection, false = "No candidate element"))


## Filter out elements overlapping their target gene of annotated promoters ------------------------

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

# overlap perturbations with features
feature_ovl <- findOverlapPairs(pert_coords, features, ignore.strand = TRUE)

# create table with overlapping features per perturbation in long format
feature_ovl <- tibble(PerturbationTargetID = first(feature_ovl)$PerturbationTargetID,
                      feature_id = second(feature_ovl)$feature_id,
                      type = second(feature_ovl)$type)

# add tested genes for perturbations overlapping any features
feature_ovl <- feature_ovl %>% 
  distinct() %>% 
  left_join(select(dl, PerturbationTargetID, measuredEnsemblID), by = "PerturbationTargetID",
            relationship = "many-to-many")

# summarize feature overlaps
feature_ovl_summary <- feature_ovl %>% 
  group_by(PerturbationTargetID, measuredEnsemblID, type) %>% 
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

# add selected overlaps to results
dl <- feature_ovl_summary %>% 
  select(PerturbationTargetID, measuredEnsemblID, ovl_target_exon, ovl_target_intron, ovl_promoter) %>% 
  left_join(dl, ., by = c("PerturbationTargetID", "measuredEnsemblID"))

# create ValidConnection column based on overlaps
dl <- dl %>% 
  mutate(ValidConnection = case_when(
    ovl_target_exon == TRUE & ValidConnection == "TRUE" ~ "overlaps target gene exon",
    ovl_target_intron == TRUE & ValidConnection == "TRUE" ~ "overlaps target gene intron",
    ovl_promoter == TRUE & ValidConnection == "TRUE" ~ "overlaps potential promoter",
    TRUE ~ ValidConnection
  ))

## Create final output -----------------------------------------------------------------------------

# manually filter out one MYC enhancer that has previously been annotated as overlapping a promoter
# and one HBG2 enhancer that isn't filtered correctly
dl <- dl %>% 
  mutate(ValidConnection = if_else(name == "MYC|chr8:128044769-128045569:.",
                                   true = "overlaps potential promoter",
                                   false = ValidConnection)) %>% 
  mutate(ValidConnection = if_else(name == "HBG2|chr11:5253147-5253547:.",
                                   true = "overlaps target gene exon",
                                   false = ValidConnection))


# reformat to ENCODE style (last 3 columns are not ENCODE format and need to be stripped for upload)
output <- dl %>% 
  select(chrom, chromStart, chromEnd, name, EffectSize, strandPerturbationTarget,
         PerturbationTargetID, chrTSS, startTSS, endTSS, strandGene,
         EffectSize95ConfidenceIntervalLow, EffectSize95ConfidenceIntervalHigh,
         measuredGeneSymbol, measuredEnsemblID, guideSpacerSeq, guideSeq, Significant, pValue,
         pValueAdjusted, starts_with("PowerAtEffectSize"), ValidConnection, Notes,
         Reference, distToTSS)

# make sure output is sorted according to genomic coordinates
output <- arrange(output, chrom, chromStart, chromEnd, measuredGeneSymbol)

# write to output
write_tsv(output, file = snakemake@output[[1]])

## Create a SingleCellExperiment object from Gasperini et al. data including perturbation status as
## alternative experiments (altExp)


# attach required packages
library(monocle)
library(SingleCellExperiment)
library(rtracklayer)
library(readr)
library(dplyr)

# load gene expression data
cds <- readRDS(snakemake@input$sce)

# load gene annotations
annot <- import(snakemake@input$annot)

# column types in guide targets file
grna_gene_pairs_cols <- cols(
  .default = col_character(),
  start.targetgene = col_double(),
  stop.targetgene = col_double()
)

# load Gasperini gRNA - gene pairs
grna_gene_pairs <- read_tsv(snakemake@input$guide_targets, col_types = grna_gene_pairs_cols)

# create SCE object from UMI counts
sce <- SingleCellExperiment(assays = list(counts = exprs(cds)))


# number of UMIs and detected genes per cell -------------------------------------------------------

# compute total umis and number of detected genes per cell
umis_per_cell <- colSums(assay(sce, "counts"))
genes_per_cell <- colSums(assay(sce, "counts") > 0)

# create DataFrame
genes_per_cell <- genes_per_cell[names(umis_per_cell)]  # anxiety check...
cell_stats <- DataFrame(total_umis = umis_per_cell,
                        detected_genes = genes_per_cell,
                        row.names = names(umis_per_cell) )

# add as colData to sce
colData(sce) <- cell_stats

# extract gene TSS annotations ---------------------------------------------------------------------

# get gene body annotations
genes <- annot[annot$type == "gene"]

# set names to gene ids without version suffix
names(genes) <- sub("\\..+", "", genes$gene_id)

## FIX THIS: only retain genes that are found in annotations
shared_genes <- intersect(rownames(sce), names(genes))
sce <- sce[shared_genes, ]
genes <- genes[shared_genes]

# get TSS coordinates
gene_tss <- resize(genes, width = 1, fix = "start")

# add TSS coordinates to SCE object as rowRanges
rowRanges(sce) <- gene_tss


# perturbation coordinates -------------------------------------------------------------------------

# remove non-targeting controls from gRNA - gene pairs
grna_gene_pairs <- grna_gene_pairs[grna_gene_pairs$general_group != "NTC", ]

# only retain data on gRNAs and their targets
grna_targets <- grna_gene_pairs %>% 
  select(chr = gRNAgroup.chr, start = gRNAgroup.start, end = gRNAgroup.stop,
         gRNAgroup, general_group) %>%
  distinct() %>% 
  mutate(start = as.integer(start), end = as.integer(end))  # convert coordinates to integers

# infer target site from gRNAgroup
grna_targets <- grna_targets %>% 
  mutate(target_site = if_else(grepl(gRNAgroup, pattern = "^chr.+\\..*"),
                               true = sub("(chr.+\\..+)_.+_.+", "\\1", gRNAgroup),
                               false = gRNAgroup))
                               
# convert to DataFrame
grna_targets <- as(grna_targets, "DataFrame")
rownames(grna_targets) <- grna_targets$gRNAgroup


# perturbation status ------------------------------------------------------------------------------

# extract cell metadata from cds and remove cds to free memory
metadat <- phenoData(cds)
rm(cds)
invisible(gc())

# extract enhancer perturbations from cell metadata
enh_perts <- varLabels(metadat)[grepl(varLabels(metadat), pattern = "^chr")]

## FIX THIS: only retain perturbations also in gRNA gene pairs
enh_perts <- intersect(enh_perts, rownames(grna_targets))

# create perturbation matrix
pert_status <- t(as(metadat[, enh_perts], "data.frame"))
pert_status <- pert_status * 1  # convert from TRUE/FALSE to 1/0
pert_status <- as(pert_status, "sparseMatrix")


# create and save output ---------------------------------------------------------------------------

# create summarized experiment with perturbation status and annotations
perts_se <- SummarizedExperiment(assays = list(perts = pert_status),
                                 rowData = grna_targets[rownames(pert_status), ])

# add perturbations to gene expression SCE object as alternative experiment
altExp(sce, e = "cre_perts") <- perts_se

# save SCE to .rds file
saveRDS(sce, file = snakemake@output[[1]])

## Create a SingleCellExperiment object from Gasperini et al. data including perturbation status as
## alternative experiments (altExp)

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages and functions
library(monocle)
library(SingleCellExperiment)
library(rtracklayer)
library(readr)
library(dplyr)
library(tidyr)
source(file.path(snakemake@scriptdir, "../R_functions/create_perturb_sce.R"))

# Load data ----------------------------------------------------------------------------------------

message("\n\nLoading data...")

# load gene expression data
cds <- readRDS(snakemake@input$cds)

# load gene annotations
annot <- import(snakemake@input$annot)

# column types in gRNA targets file
grna_targets_cols <- cols(
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

# load gRNA targets
grna_targets <- read_tsv(snakemake@input$guide_targets, col_types = grna_targets_cols)

# extract gene expression and column data from cds object
expr <- exprs(cds)
col_data <- phenoData(cds)

# remove cds to free memory
rm(cds)
invisible(gc())

# get gene TSS annotations -------------------------------------------------------------------------

message("Getting TSS annotations...")

# extract gene locus annotations
genes <- annot[annot$type == "gene"]

# set names to gene ids without version suffix
names(genes) <- sub("\\..+", "", genes$gene_id)

# get TSS coordinates
gene_tss <- resize(genes, width = 1, fix = "start")

# check if any genes in expression data are missing from gene annotations
missing_genes <- length(setdiff(rownames(expr), names(genes)))
if (length(missing_genes) > 0) {
  warning(missing_genes, " genes in expr not found in annot will be dropped.", call. = FALSE)
}

# only retain expression for genes that are also found in annotations
expr <- expr[intersect(rownames(expr), names(genes)), ]

# extract detected gRNAs per cell ------------------------------------------------------------------

message("Getting detected gRNAs per cell...")

# extract gRNAs per cell
grnas_per_cell <- as(col_data[, c("cell", "barcode")], "data.frame")

# split barcode string in to individual guides and convert to long format
grnas_per_cell <- grnas_per_cell %>% 
  as_tibble() %>% 
  mutate(barcode = strsplit(barcode, split = "_")) %>% 
  unnest(cols = barcode) %>% 
  rename(spacer = barcode)

# add gRNA id and remove sequence, since that is not required
grnas_per_cell <- grnas_per_cell %>% 
  left_join(distinct(select(grna_targets, name, spacer)), by = "spacer") %>% 
  select(-spacer) %>% 
  filter(!is.na(name)) %>% 
  rename(grna = name)

# create output ------------------------------------------------------------------------------------

message("Creating SCE object...")

# create Perturb-seq SingleCellExperiment object
sce <- create_pert_sce(expr, annot = gene_tss, pert_status = grnas_per_cell,
                       grna_annot = grna_targets, targets_pert_name = "cre_perts")

# save SCE to .rds file
message("Writing SCE object to file...")
saveRDS(sce, file = snakemake@output[[1]])

message("Done!")

# close log file connection
sink()
sink(type = "message")
close(log)

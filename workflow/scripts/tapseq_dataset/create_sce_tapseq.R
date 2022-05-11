## Create a Perturb-seq SingleCellExperiment object from TAP-seq processing workflow output

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(rtracklayer)
  library(readr)
  source(file.path(snakemake@scriptdir, "../R_functions/create_perturb_sce.R"))
})

# column types in guide targets file
guide_targets_cols <- cols(
  .default = col_character(),  # for target_type columns, which is optional as of now
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  strand = col_character(),
  spacer = col_character(),
  target_chr = col_character(),
  target_start = col_integer(),
  target_end = col_integer(),
  target_name = col_character(),
  target_strand = col_character()
)

# load guide targets
guide_targets <- read_tsv(snakemake@input$guide_targets, col_types = guide_targets_cols)

# process perturbation status matrix ---------------------------------------------------------------

message("Processing perturbation status matrix...")

# load perturbation status matrix
pert_status <- fread(snakemake@input$pert_status)

# convert to long format
pert_status <- melt(pert_status, id.vars = "VECTOR", variable.name = "cell", value.name = "pert",
                    variable.factor = FALSE)

# only retain observed cell-guide combinations as required by the create_pert_sce() function
pert_status <- pert_status[pert_status$pert > 0, c("cell", "VECTOR")]
colnames(pert_status) <- c("cell", "grna")

# filter for guides in guide targets and convert to data.frame for create_pert_sce() function
pert_status <- as.data.frame(pert_status[pert_status$grna %in% guide_targets$name, ])

# free up memory
message("Freeing up memory...")
gc()

# process dge matrix -------------------------------------------------------------------------------

message("Processing DGE matrix...")

# load dge matrix
expr <- fread(snakemake@input$dge)

# remove any data on vector expression
vectors <- grepl(expr[[1]], pattern = snakemake@params$vector_pattern)
expr <- expr[!vectors, ]

# convert to sparse matrix
expr <- as(as.matrix(expr, rownames = 1), "sparseMatrix")

# get gene TSS annotations -------------------------------------------------------------------------

message("Getting TSS annotations...")

# load gene annotations
annot <- import(snakemake@input$annot)

# extract gene locus annotations and set names to gene names
genes <- annot[annot$type == "gene"]
names(genes) <- genes$gene_name

# calculate TSS coordinates
gene_tss <- resize(genes, width = 1, fix = "start")

# check if any genes in expression data are missing from gene annotations
missing_genes <- length(setdiff(rownames(expr), names(genes)))
if (missing_genes > 0) {
  warning(missing_genes, " genes in expr not found in annot will be dropped.", call. = FALSE)
}

# only retain expression for genes that are also found in annotations
expr <- expr[intersect(rownames(expr), names(genes)), ]

# create Perturb-seq SCE object --------------------------------------------------------------------

message("Creating SCE object...")

# create Perturb-seq SingleCellExperiment object
sce <- create_pert_sce(expr, annot = gene_tss, pert_status = pert_status,
                       grna_annot = guide_targets, targets_pert_name = "cre_perts")

# save SCE to .rds file
message("Writing SCE object to file...")
saveRDS(sce, file = snakemake@output[[1]])

message("Done!")

# close log file connection
sink()
sink(type = "message")
close(log)

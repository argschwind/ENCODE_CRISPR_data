## This snakemake R-script performs differential gene expression tests for trans-effects in
## Perturb-seq data

# save.image(paste0("trans_effects.rda"))
# stop()

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages and functions
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(BiocParallel)
  library(SingleCellExperiment)
  source(file.path(snakemake@scriptdir, "../R_functions/differential_expression_fun.R"))
  source(file.path(snakemake@scriptdir, "trans_differential_expression_fun.R"))
})

# parse method wildcard and attach required packages
method <- snakemake@wildcards$method
library(method, character.only = TRUE)

# register parallel backend if specified (if more than 1 thread provided) and set RNG seed
if (snakemake@threads > 1) {
  message("Registering parallel backend with ", snakemake@threads, " cores.")
  register(MulticoreParam(workers = snakemake@threads,
                          RNGseed = snakemake@params$seed))
} else {
  message("Registering serial backend.")
  register(SerialParam(RNGseed = snakemake@params$seed))
}

# prepare data =====================================================================================

# load prepared input data stored in SingleCellExperiment object
message("Loading input data.")
sce <- readRDS(snakemake@input$sce)

# infer perturbation level based on strategy
pert_level <- switch(snakemake@wildcards$strategy, "perGRNA" = "grna_perts", "perCRE" = "cre_perts",
                     stop("incorrect strategy argument"))

# filter cells for minimum and maximum number total UMIs per cell
sce <- filter_umis_per_cell(sce, min_umis = snakemake@params$umis_per_cell[[1]],
                            max_umis = snakemake@params$umis_per_cell[[2]])

# filter for minimum number of cells per perturbation
sce <- filter_cells_per_pert(sce, min_cells = snakemake@params$min_cells, pert_level = pert_level)

# perform DE tests =================================================================================

# normalize counts using censored mean normalization and log transformation
message("Normalizing transcript counts.")
invisible(gc())
sce <- normalize_cens_mean(sce)
invisible(gc())
assay(sce, "logcounts") <- log1p(assay(sce, "normcounts"))
invisible(gc())

# get differential expression function based on method from wildcards
de_function <- get(paste0("de_", method))

# get perturbation to analyze in this batch
pert_batches <- read_tsv(snakemake@input$pert_batches, show_col_types = FALSE)
perts <- pull(filter(pert_batches, batch == snakemake@wildcards$batch), pert)

# only retain perturbations that were not filtered out due to the low cell numbers
filtered_perts <- rownames(assay(altExp(sce, pert_level), "perts"))
perts <- intersect(perts, filtered_perts)
names(perts) <- perts

# perform differential gene expression analysis
message("Performing differential expression tests.")
output <- test_trans_diff_expression(sce, perts = perts,
                                     pert_level = pert_level,
                                     sample_genes = snakemake@params$sample_genes,
                                     de_function = de_function,
                                     formula = snakemake@params$formula,
                                     n_ctrl = snakemake@params$n_ctrl,
                                     cell_batches = snakemake@params$cell_batches,
                                     p_adj_method = snakemake@params$p_adj_method)

# reformat output
message("Reformating output.")
output <- output %>% 
  select(-c(gene_end, pert_center)) %>% 
  dplyr::rename(gene_tss = gene_start, dist_to_tss = distance)

# save DE output to file
message("Saving output to file.")
write_tsv(output, file = snakemake@output[[1]])

# close log file connection
sink()
sink(type = "message")
close(log)

## This snakemake R-script performs differential gene expression analysis for Perturb-seq
## discover regulatory interactions from perturbation - gene pairs

# save.image(paste0("diff_expr_", snakemake@wildcards$chr, ".rda"))
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
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R")) 
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
sce <- readRDS(snakemake@input[[1]])

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
sce <- normalize_cens_mean(sce)
assay(sce, "logcounts") <- log1p(assay(sce, "normcounts"))

# ensure that max distance if numeric if not NULL
if (!is.null(snakemake@params$max_dist)) {
  max_dist <- as.numeric(snakemake@params$max_dist)
} else {
  max_dist <- snakemake@params$max_dist
}

# get differential expression function based on method from wildcards
de_function <- get(paste0("de_", method))

# perform differential gene expression analysis
message("Performing differential expression tests.")
output <- test_differential_expression(sce, pert_level = pert_level,
                                       max_dist = max_dist,
                                       de_function = de_function,
                                       formula = snakemake@params$formula,
                                       n_ctrl = snakemake@params$n_ctrl,
                                       cell_batches = snakemake@params$cell_batches,
                                       p_adj_method = snakemake@params$p_adj_method)

# reformat output
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

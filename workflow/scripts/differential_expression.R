## This snakemake R-script performs differential gene expression analysis for Perturb-seq
## perturbations to discover cis-regulatory interactions


# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages and functions
suppressPackageStartupMessages({
  library(dplyr)
  library(BiocParallel)
  library(SingleCellExperiment)
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R")) 
})

# register parallel backend if specified (if more than 1 thread provided)
if (snakemake@threads > 1) {
  register(MulticoreParam(workers = snakemake@threads))
} else {
  register(SerialParam())
}

# set seed for reproducible results
if (snakemake@params$seed > 0) set.seed(snakemake@params$seed)

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

# perform differential gene expression analysis
message("Performing differential expression tests.")
output <- test_differential_expression(sce, pert_level = pert_level,
                                       max_dist = as.numeric(snakemake@params$max_dist),
                                       method = snakemake@wildcards$method,
                                       formula = snakemake@params$formula,
                                       n_ctrl = snakemake@params$n_ctrl,
                                       cell_batches = snakemake@params$cell_batches)

# reformat output
output <- output %>% 
  select(-c(gene_chr, gene_end, pert_center)) %>% 
  dplyr::rename(chr = pert_chr, gene_tss = gene_start, dist_to_tss = distance)

# save DE output to file
write.csv(output, file = snakemake@output[[1]], row.names = FALSE)

# close log file connection
sink()
sink(type = "message")
close(log)

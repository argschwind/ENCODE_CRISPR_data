## Fit negative binomial distributions for each gene using DESeq2

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages and functions
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R"))
})

# load prepared input data stored in SingleCellExperiment object
message("Loading input data.")
sce <- readRDS(snakemake@input[[1]])

# filter cells for minimum and maximum number total UMIs per cell
sce <- filter_umis_per_cell(sce, min_umis = snakemake@params$umis_per_cell[[1]],
                            max_umis = snakemake@params$umis_per_cell[[2]])

# remove specific genes from power simulations
sce <- sce[!rownames(sce) %in% snakemake@params$remove_genes, ]

# fit negative binomial distributions to estimate gene-level dispersion
message("Estimate dispersion using DESeq2:")
sce <- fit_negbinom_deseq2(sce, size_factors = snakemake@params$size_factors,
                           fit_type = snakemake@params$fit_type)

# save output sce to file
saveRDS(sce, file = snakemake@output[[1]])

# close log file connection
sink()
sink(type = "message")
close(log)

# Perform power simulations for Perturb-seq style experiments. For each perturbation a specified
# decrease in expression is simulated for each gene and differential expression tests are performed
# to assess if the expression change can be detected, i.e. if the data for a given perturbation has
# enough power to detect a change of a given effect size.

#save.image(paste0("pwr_sim_", snakemake@wildcards$chr, ".rda"))
#stop()

# opening log file to collect all messages, warnings and errors
tmp_logfile <- tools::file_path_sans_ext(snakemake@log[[1]])
log <- file(tmp_logfile, open = "w")
sink(log)
sink(log, type = "message")

# required packages and functions
suppressPackageStartupMessages({
  library(dplyr)
  library(BiocParallel)
  library(SingleCellExperiment)
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
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

# remove specific genes from power simulations
sce <- sce[!rownames(sce) %in% snakemake@params$remove_genes, ]

# perform power simulations ========================================================================

# normalize counts using censored mean normalization and log transformation
message("Normalizing transcript counts.")
sce <- normalize_cens_mean(sce)
assay(sce, "logcounts") <- log1p(assay(sce, "normcounts"))

# fit negative binomial distributions to estimate gene-level dispersion
message("Estimate dispersion using DESeq2:")
sce <- fit_negbinom_deseq2(sce, size_factors = snakemake@params$size_factors,
                           fit_type = snakemake@params$fit_type)

# convert 'percentage decrease' effect size to 'relative expression level'
effect_size <- 1 - as.numeric(snakemake@wildcards$effect)

# simulate TAP-seq data and perform differential gene expression tests
message("Performing power simulations.")
output <- simulate_diff_expr(sce, effect_size = effect_size,
                             pert_level = pert_level,
                             max_dist = as.numeric(snakemake@params$max_dist),
                             genes_iter = snakemake@params$pert_genes,
                             guide_sd = as.numeric(snakemake@wildcards$sd),
                             center = FALSE,
                             rep = snakemake@params$rep,
                             method = snakemake@wildcards$method,
                             formula = as.formula(snakemake@params$formula),
                             n_ctrl = snakemake@params$n_ctrl,
                             cell_batches = snakemake@params$cell_batches)

# extract DESeq2 outlier information for every gene
message("Processing output.")
disp_outlier <- data.frame(gene = rownames(rowData(sce)),
                           disp_outlier_deseq2 = rowData(sce)[, "disp_outlier_deseq2"],
                           stringsAsFactors = FALSE)

# add to output
output <- left_join(output, disp_outlier, by = "gene")

# save simulation output
message("Saving output to file.")
outfile <- gzfile(snakemake@output[[1]], open = "w")
write.csv(output, file = outfile, row.names = FALSE)
close(outfile)

# close and compress log file
sink()
sink(type = "message")
close(log)
R.utils::gzip(tmp_logfile)

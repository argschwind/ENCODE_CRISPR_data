## Infer 10x lane barcode from cell barcodes in TAP-seq data and add to colData of Perturb-seq
## SingleCellExperiment object. Filter out cells belonging to specified lanes.

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages and functions
suppressPackageStartupMessages({
  library(SingleCellExperiment)
}) 

# load Perturb-seq SCE object
sce <- readRDS(snakemake@input[[1]])

# get 10x lane barcodes from cell barcodes and add to colData
cell_barcodes <- colnames(sce)
lane_barcodes <- substr(cell_barcodes, start = 1, stop = 8)
colData(sce)[, "lane_barcodes"] <- lane_barcodes

# filter out cells from specified 10x lanes
lane_filter <- colData(sce)[, "lane_barcodes"] %in% snakemake@params$remove_lanes
sce <- sce[, !lane_filter]

# report number of filtered cells
message("Filtering out ", sum(lane_filter), " cells from the following 10x lanes:\n",
        paste(snakemake@params$remove_lanes, collapse = ", "))

# write filtered SCE object to output
saveRDS(sce, file = snakemake@output[[1]])

# close log file connection
sink()
sink(type = "message")
close(log)

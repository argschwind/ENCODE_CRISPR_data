## Extract data on one chromosome from a Perturb-seq (or TAP-seq) SingleCellExperiment object and
## save to new file. Used to split data into junks by chromosome to parallelize mapping regulatory
## interactions

# attach required packages
suppressPackageStartupMessages(library(SingleCellExperiment))

# load input SCE
sce <- readRDS(snakemake@input[[1]])

# extract data for genes on selected chromosome
chr <- snakemake@wildcards$chr
sce <- sce[seqnames(rowRanges(sce)) == chr, ]

# helper function to subset an alternative experiment based on chromosome. requires that rowData has
# a column called "chr" with the chromosome id of that feature
subset_altexp <- function(sce, chr) {
  
  # all alternative experiments in sce
  alt_exps <- altExpNames(sce)
  
  # subset all altExps in SCE by provided chromosome ids
  for (e in alt_exps) {
    alt <- altExp(sce, e = e)
    alt_chr <- alt[rowData(alt)[["chr"]] %in% chr, ]
    altExp(sce, e = e) <- alt_chr
  }
  
  return(sce)
  
}

# subset perturbation data in sce based on chromosome
sce <- subset_altexp(sce, chr = chr)

# remove any cells with 0 zero counts for this chromosome if specified
if (snakemake@params$rm_zero_cells == TRUE) {
  
  # get cells with 0 counts for genes on chromosome
  zeros <- colSums(assay(sce, "counts")) == 0
  
  # remove these cells from sce
  message("\nFiltering out ", sum(zeros), " cells with 0 counts for genes on chromosome ", chr, ".")
  sce <- sce[, !zeros]
  
}

# save to output .rds file
saveRDS(sce, file = snakemake@output[[1]])

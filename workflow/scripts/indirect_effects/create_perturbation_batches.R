
# required packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
})

# load SCE object
sce <- readRDS(snakemake@input[[1]])

# infer perturbation level based on strategy
pert_level <- switch(snakemake@wildcards$strategy, "perGRNA" = "grna_perts", "perCRE" = "cre_perts",
                     stop("incorrect strategy argument"))

# extract all perturbations
perts <- rownames(altExp(sce, pert_level))

# create table with all perturbations and randomly assign them to 1 batch
n_batches <- snakemake@params$n_batches
batches <- rep(1:n_batches, times = ceiling(length(perts) / n_batches))
pert_batches <- perts %>% 
  tibble(pert = .) %>% 
  mutate(batch = batches[1:n()])

# write batches to output file
write_tsv(pert_batches, file = snakemake@output[[1]])

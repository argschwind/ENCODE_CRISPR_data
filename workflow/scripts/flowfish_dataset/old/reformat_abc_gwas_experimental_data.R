## Reformat Fulco and Gasperini FlowFish data for ENCODE4: Lift enhancer coordinates from hg19 to
## hg38, reformat into ENCODE E-G benchmarking experimental data format and split into Fulco and 
## Fulco & Gasperini data sets.

library(tidyverse)
library(here)
library(rtracklayer)


## Load data ---------------------------------------------------------------------------------------

# hg19 to hg38 file used for liftover performed by this script
chain_url <- "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

# download chain file
chain_file <- file.path(tempdir(), basename(chain_url))
download.file(chain_url, chain_file)
system(paste("gunzip", chain_file))  # de-compress, since rtracklayer can't import compressed chains

# import liftover chain
hg19_to_hg38_chain <- import.chain(tools::file_path_sans_ext(chain_file))

# url to file containing all experimental data from GWAS ABC paper
expt_url <- "https://raw.githubusercontent.com/EngreitzLab/ABC-GWAS-Paper/main/comparePredictorsToCRISPRData/comparisonRuns/K562-only/experimentalData/experimentalData.K562-only.txt"

# column types in expt file
expt_cols <- cols(
  .default = col_character(),
  startPerturbationTarget = col_double(),
  endPerturbationTarget = col_double(),
  startTSS = col_double(),
  endTSS = col_double(),
  EffectSize = col_double(),
  Significant = col_logical(),
  padj = col_double(),
  Regulated = col_logical(),
  IncludeInModel = col_logical(),
  PowerAtEffectSize.0.25 = col_double()
)

# load experimental data
expt <- read_tsv(expt_url, col_types = expt_cols)


## Reformat data to ENCODE format ------------------------------------------------------------------

# rename power columns
power_cols <- grepl(colnames(expt), pattern = "PowerAtEffectSize")
effect_sizes <- as.numeric(sub("PowerAtEffectSize.", "", colnames(expt)[power_cols])) * 100
colnames(expt)[power_cols] <- paste0("PowerAtEffectSize", effect_sizes)

# create ValidConnection column
expt <- expt %>% 
  mutate(ValidConnection = if_else(IncludeInModel == TRUE, true = "TRUE", false = Reason))

# create name column based on enhancer coordinates and target gene
expt <- mutate(expt, name = paste0(GeneSymbol, "|", chrPerturbationTarget, ":",
                                   startPerturbationTarget, "-", endPerturbationTarget, ":*"))

# extract required columns for hg19 experimental data and sort according to coordinates
expt_hg19 <- expt %>% 
  select(chrom = chrPerturbationTarget, chromStart = startPerturbationTarget,
         chromEnd = endPerturbationTarget, name, EffectSize, chrTSS, startTSS, endTSS,
         measuredGeneSymbol = GeneSymbol, Significant, pValueAdjusted = padj,
         PowerAtEffectSize25, ValidConnection, CellType, Reference, Regulated) %>% 
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)


## Liftover enhancer coordinates -------------------------------------------------------------------

# add unique identifier for enhancers
expt_hg19 <- mutate(expt_hg19, enh_uid = paste(chrom, chromStart, chromEnd, sep = "_"))

# extract hg19 enhancer coordinates and create GenomicRanges object
enh_coords_hg19 <- expt_hg19 %>%
  select(chrom, chromStart, chromEnd, enh_uid) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# liftover enhancers from hg19 to hg38
enh_coords_hg38 <- liftOver(enh_coords_hg19, chain = hg19_to_hg38_chain)

# get cases where one enhancer was mapped to multiple hg38 locations
n_hg38 <- vapply(enh_coords_hg38, FUN = length, FUN.VALUE = integer(1))
multiple_hg38 <- enh_coords_hg38[n_hg38 > 1] %>% 
  unlist() %>% 
  mcols() %>%
  as_tibble() %>% 
  pull(enh_uid) %>% 
  unique()

# get experimental data involving enhancers with hg38 multi-mapping issues and save to file
expt_hg38_multimap <- filter(expt_hg19, enh_uid %in% multiple_hg38)
mm_outfile <- paste0(tools::file_path_sans_ext(basename(expt_url)), "_multi_hg38_mappings.tsv")
write_tsv(select(expt_hg38_multimap, -enh_uid),
          file = here("reformatted_experimental_data", mm_outfile))

# set enh_uids as names for hg38 enhancer coordinates GRangesList
names(enh_coords_hg38) <- vapply(enh_coords_hg38, FUN = function(x) {
  unique(mcols(x)[["enh_uid"]])
}, FUN.VALUE = character(1))

# function to resolve cases where enhancers mapped to multiple hg38 locations by selecting the
# mapping where the width of the mapped location is closest to the original width
resolve_multi_mapping <- function(mapped, input) {
  mapped <- mapped[seqnames(mapped) == seqnames(input)]  # prevent liftover to other chromosomes
  mapped[which.min(abs(width(input) -  width(mapped)))]
}

# apply function
enh_coords_hg19 <- split(enh_coords_hg19, f = enh_coords_hg19$enh_uid)
enh_coords_hg38 <- enh_coords_hg38[names(enh_coords_hg19)]
enh_coords_hg38_1mapped <- mapply(FUN = resolve_multi_mapping, enh_coords_hg38, enh_coords_hg19)

# convert hg38 coordinates to data frame
hg38_coords_df <- enh_coords_hg38_1mapped %>% 
  GRangesList() %>% 
  unlist() %>% 
  as_tibble() %>% 
  select(chrom = seqnames, chromStart = start, chromEnd = end, enh_uid) %>% 
  mutate(chrom = as.character(chrom))

# replace hg19 coordinates with hg38 coordinates
expt_cols <- colnames(expt_hg19) 
expt_hg38 <- expt_hg19 %>%
  select(-c(chrom, chromStart, chromEnd)) %>%
  left_join(hg38_coords_df, by = "enh_uid") %>% 
  select(!!expt_cols)


## Create output files -----------------------------------------------------------------------------

# remove enhancer uid column
expt_hg19 <- select(expt_hg19, -enh_uid)
expt_hg38 <- select(expt_hg38, -enh_uid)

# filter based on 80% power to detect 25% expression change
expt_hg19_filt <- filter(expt_hg19, Significant == TRUE | PowerAtEffectSize25 >= 0.8)
expt_hg38_filt <- filter(expt_hg38, Significant == TRUE | PowerAtEffectSize25 >= 0.8)

# write to output files
outdir <- "reformatted_experimental_data"
write_tsv(expt_hg19_filt, file = here(outdir, "ExperimentalData_ABC-GWAS_K562_hg19_210712.tsv.gz"))
write_tsv(expt_hg38_filt, file = here(outdir, "ExperimentalData_ABC-GWAS_K562_GRCh38_210712.tsv.gz"))

## Liftover CRISPR enhancer screen datasets from hg19 to hg38

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rtracklayer)
})

# Load data ----------------------------------------------------------------------------------------

# relevant column types in CRISPR data (can be either full ENCODE data or EP benchmarking subset)
crispr_data_cols <- cols(
  chrom = col_character(),
  chromStart = col_integer(),
  chromEnd = col_integer(),
  name = col_character(),
  chrTSS = col_character(),
  startTSS = col_integer(),
  endTSS = col_integer()
)

# load CRISPR dataset in hg19
dat <- read_tsv(snakemake@input$results, col_types = crispr_data_cols)

# column types in a 6 column bed file (score as character due to '.' for empty )
bed_cols <- cols(
  chrom = col_character(),
  chromStart = col_integer(),
  chromEnd = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character()
)

# load enhancer coordinates lifted over to hg38
enh_hg38 <- read_tsv(snakemake@input$enh_hg38, col_names = names(bed_cols$cols),
                     col_types = bed_cols)

# load hg38 annotations
annot_hg38 <- import(snakemake@input$annot_hg38, format = "gtf")

# Liftover enhancer coordinates --------------------------------------------------------------------

# get original column names of dat
dat_colnames <- colnames(dat)

# add unique identifier for enhancers
dat <- mutate(dat, enh_uid = paste0(chrom, ":", chromStart, "-", chromEnd))

# replace hg19 enhancer coordinates by hg38 coordinates
dat <- dat %>% 
  select(-c(chrom, chromStart, chromEnd)) %>% 
  left_join(select(enh_hg38, -c(score, strand)), by = c("enh_uid" = "name"))

# remove any pairs involving enhancers that failed to lift over
unlifted_filter <- is.na(dat$chrom)
if (any(unlifted_filter)) {
  warning("Removing ", sum(unlifted_filter), " E-G pairs for which enhancer liftover failed!",
          call. = FALSE)
  dat <- dat[!unlifted_filter, ]
}

# Liftover TSS coordinates -------------------------------------------------------------------------

# only retain annotations on genes in data
annot_hg38 <- annot_hg38[annot_hg38$gene_name %in% unique(dat$measuredGeneSymbol)]

# extract gene locus coordinates
genes_hg38 <- annot_hg38[annot_hg38$type == "gene"]

# merge weird gene loci with multiple annotations
genes_hg38 <- split(genes_hg38, f = genes_hg38$gene_name)
genes_hg38 <- reduce(genes_hg38)
genes_hg38 <- unlist(genes_hg38)

# create TSS annotations based on locus coordinates (0-based)
tss_hg38 <- resize(genes_hg38, width = 2, fix = "start")

# convert to data frame with data to be added to hg19 data
tss_hg38 <- data.frame(chrTSS = as.character(seqnames(tss_hg38)),
                       startTSS = start(tss_hg38),
                       endTSS = end(tss_hg38),
                       measuredGeneSymbol = names(tss_hg38),
                       stringsAsFactors = FALSE)

# replace hg19 TSS coordinates in data with hg38 coordinates
output <- dat %>% 
  select(-c(startTSS, endTSS)) %>% 
  left_join(tss_hg38, by = c("measuredGeneSymbol", "chrTSS"))

# rearrange columns to original order for output
output <- output %>% 
  select(all_of(dat_colnames)) %>% 
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)
  
# write output to file
write_tsv(output, file = snakemake@output[[1]])

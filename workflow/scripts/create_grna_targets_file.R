## Create a guide targets file based on gRNA binding coordinates and candidate target coordinates.
## Overlaps gRNAs with targets and assign each guide to overlapping targets.

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
library(dplyr)
library(readr)
library(rtracklayer)

# Load data ----------------------------------------------------------------------------------------

message("\n\nLoading data...")

# types of the first 7 columns in a gRNAs bed file
grnas_cols <- cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character(),
  spacer = col_character()
)

# column names of the first 7 columns in a bed file
grnas_col_names <- c("chr", "start", "end", "name", "score", "strand", "spacer")

# load gRNA positions
grnas <- read_tsv(snakemake@input$grnas, col_names = grnas_col_names, col_types = grnas_cols)
grnas <- grnas[, c(1:4, 6, 7)]  # only retain needed columns

# types of the first 6 columns in a standard bed file
bed_cols <- cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character(),
)

# column of the first 6 columns in a standard bed file
bed_col_names <- grnas_col_names[1:6]

# load candidate target elements coordinates
targets <- read_tsv(snakemake@input$targets, col_names = bed_col_names, col_types = bed_cols)
targets <- targets[, c(1:4, 6)]  # only retain needed columns

# Create guide targets file ------------------------------------------------------------------------

message("Creating guide targets file...")

# create GRanges objects
grnas_gr <- makeGRangesFromDataFrame(grnas, keep.extra.columns = TRUE,
                                     starts.in.df.are.0based = TRUE)
targets_gr <- makeGRangesFromDataFrame(targets, keep.extra.columns = TRUE,
                                       starts.in.df.are.0based = TRUE)

# check if target regions overlap
merged_targets <- reduce(targets_gr)
if (length(merged_targets) < length(targets_gr)) {
  warning("Some targets overlap. gRNAs targeting those might get assiged to multiple targets!",
          call. = FALSE)
}

# find overlaps
ovl <- findOverlaps(query = grnas_gr, subject = targets_gr)

# extract gRNAs and targets of overlapping pairs
grna_targets <- cbind(grnas[queryHits(ovl), ], targets[subjectHits(ovl), ])
colnames(grna_targets)[7:11] <- paste0("target_", colnames(grna_targets)[7:11])

# count the number of overlaps for every gRNA, including those not overlapping any region
ovl_counts <- grna_targets %>% 
  mutate(name = factor(name, levels = unique(grnas$name))) %>% 
  count(name, .drop = FALSE)

# save output to files
write_tsv(grna_targets, file = snakemake@output$grna_targets_file)
write_tsv(ovl_counts, file = snakemake@output$n_overlaps)

message("Done!")

# close log file connection
sink()
sink(type = "message")
close(log)

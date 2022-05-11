## Create a guide targets file based on gRNA binding positions and candidate target coordinates.
## Overlaps guides with targets and assign each guide to overlapping targets. Includes any known
## guide targets specified via a known_targets file

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rtracklayer)
})

# Load data ----------------------------------------------------------------------------------------

message("\n\nLoading data...")

# types of the first 7 columns in a guide bed file
guide_cols <- cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character(),
  spacer = col_character()
)

# load guide positions and only retain needed columns
guides <- read_tsv(snakemake@input$guides, col_types = guide_cols, col_names = names(guide_cols$cols))
guides <- guides[, c(1:4, 6, 7)]
guide_ids <- guides$name  # extract all guide ids for later use

# types of the first 6 columns in a standard bed file
bed_cols <- cols(
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  score = col_double(),
  strand = col_character(),
)

# load candidate target elements coordinates and only retain needed columns
targets <- read_tsv(snakemake@input$targets, col_types = bed_cols, col_names = names(bed_cols$cols))
targets <- targets[, c(1:4, 6)]

# load known targets file if specified and extract required columns
if (!is.null(snakemake@input$known_targets)) {
  
  # load known targets files
  known_targets <- read_tsv(snakemake@input$known_targets, show_col_types = FALSE)
  
  # extract required columns in correct order (any other columns are ignored)
  known_target_cols <- c("name", "target_chr", "target_start", "target_end", "target_name",
                         "target_strand", "target_type")
  known_targets <- select(known_targets, all_of(known_target_cols))
  
}

# Create guide targets file ------------------------------------------------------------------------

message("Creating guide targets file...")

# add guide positions to known targets and remove these from guides, since they do not need to be
# overlapped with the provided target elements
if (!is.null(snakemake@input$known_targets)) {
  
  message("Processing known targets...")
  
  # add guide positions to known targets
  known_targets <- right_join(guides, known_targets, by = "name")
  
  # remove known targets guides from guide list
  guides <- filter(guides, !name %in% known_targets$name)
  
}

# create GRanges objects
guides_gr <- makeGRangesFromDataFrame(guides, starts.in.df.are.0based = TRUE)
targets_gr <- makeGRangesFromDataFrame(targets, starts.in.df.are.0based = TRUE)

# check if target regions overlap
merged_targets <- reduce(targets_gr)
if (length(merged_targets) < length(targets_gr)) {
  warning("Some targets overlap. gRNAs targeting those might get assiged to multiple targets!",
          call. = FALSE)
}

# find overlaps
ovl <- findOverlaps(query = guides_gr, subject = targets_gr, ignore.strand = TRUE)

# extract guides and targets of overlapping pairs
colnames(targets) <- paste0("target_", colnames(targets))
guide_targets <- bind_cols(guides[queryHits(ovl), ], targets[subjectHits(ovl), ])

# add target_type for guides that were overlapped with candidate targets based on provided type
guide_targets <- mutate(guide_targets, target_type = snakemake@params$targets_type)

# add known targets and sort the same as input guide file
guide_targets <- bind_rows(known_targets, guide_targets) %>% 
  mutate(name = factor(name, levels = guide_ids)) %>% 
  arrange(name)

# count the number of targets per guide, including those not overlapping any target element
targets_per_guide <- count(guide_targets, name, name = "targets", .drop = FALSE)

# create data frame in bed format with target coordinates
targets_bed <- guide_targets %>% 
  mutate(score = if_else(target_type == "known", true = 1000, false = 500)) %>% 
  select(target_chr, target_start, target_end, target_name, score, target_strand) %>% 
  distinct() %>% 
  arrange(target_chr, target_start, target_end, target_name)

# save output to files
write_tsv(guide_targets, file = snakemake@output$guide_targets_file)
write_tsv(targets_per_guide, file = snakemake@output$targets_per_guide_file)
write_tsv(targets_bed, file = snakemake@output$targets_bed, col_names = FALSE)

message("Done!")

# close log file connection
sink()
sink(type = "message")
close(log)

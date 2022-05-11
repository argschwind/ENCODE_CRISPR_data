## Match TAP-seq guides to hg19 DNA sequence of target (if available) to identify binding sites of
## guides

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(Biostrings)
  library(BSgenome)
  library(BiocParallel)
})

# register backend for parallel computing
if (snakemake@threads > 1) {
  register(MulticoreParam(workers = snakemake@threads))
} else {
  register(SerialParam())
}

# define functions =================================================================================

# get guide sequence from vector sequence
get_guide_seq <- function(vector, upstream_seq, downstream_seq) {
  
  # get end of upstream sequence and start of downstream sequence
  end_upstream <- end(matchPattern(vector, pattern = upstream_seq))
  start_downstream <- start(matchPattern(vector, pattern = downstream_seq))
  
  # extract guide sequence
  guide_seq <- subseq(vector, start = end_upstream + 1, end = start_downstream - 1)
  
  return(guide_seq)
  
}

# function to get coordinates for one guide binding site
get_guide_coords <- function(guide_id, guide_seqs, guide_chrs, genome, pam = "NGG") {
  
  # get guide sequence and chromosome
  seq <- guide_seqs[[guide_id]]
  chr <- guide_chrs[guide_chrs$guide == guide_id, "chr"]
  
  # add PAM site to guide sequence
  seq_pam <- xscat(seq, pam)
  
  # get binding sites for positive and negative strand
  match_pos <- matchPattern(seq_pam, subject = genome[[chr]], fixed = "subject")
  match_neg <- matchPattern(reverseComplement(seq_pam), subject = genome[[chr]], fixed = "subject")
  
  # create GRanges with coordinates of binding site on positive strand
  if (length(match_pos)) {
    guide_ranges <- ranges(match_pos)
    end(guide_ranges) <- end(guide_ranges) - nchar(pam)
    coords_pos <- GRanges(seqnames = chr, ranges = guide_ranges, strand = "+",
                          name = rep(guide_id, length(match_pos)),
                          seq = rep(as.character(seq), length(match_pos)))
  } else {
    coords_pos <- GRanges()
  }
  
  # create GRanges with coordinates of binding site on negative strand
  if (length(match_neg)) {
    guide_ranges <- ranges(match_neg)
    start(guide_ranges) <- start(guide_ranges) + nchar(pam)
    coords_neg <- GRanges(seqnames = chr, ranges = guide_ranges, strand = "-",
                          name = rep(guide_id, length(match_neg)),
                          seq = rep(as.character(seq), length(match_neg)))
  } else {
    coords_neg <- GRanges()
  }
  
  # combine to create output
  guide_coords <- c(coords_pos, coords_neg)
  
  return(guide_coords)
  
}

# retain only guide binding sites within a given target site
filter_guide_coords_target <- function(guide_coord, guide_targets) {
  target_sites <- subsetByOverlaps(guide_coord, guide_targets[[unique(guide_coord$name)]])
  if (length(target_sites) > 1) {
    warning("Multiple binding site within target found for: ", unique(guide_coord$name),
            ". Returning only first site.", call. = FALSE)
  }
  return(target_sites[1])
}

# create targets for control guides
ctrl_guide_targets <- function(guides) {
  target <- rep(range(guides, ignore.strand = TRUE), times = length(guides))
  target$guide_id <- guides$name
  return(target)
}

# create guide and target annotations ==============================================================

# load vector sequences
vectors <- readDNAStringSet(snakemake@input$vectors, format = "fasta")

# load target gene annotations and extract gene annotations
annot <- import(snakemake@input$annot, format = "gtf")
genes <- annot[annot$type == "gene"]

# get chromosome for each gene
gene_chrs <- data.frame(gene = genes$gene_name, chr = as.character(seqnames(genes)),
                        stringsAsFactors = FALSE)
gene_chrs <- distinct(gene_chrs)

# get guide sequences ------------------------------------------------------------------------------

# get guide sequences from all vectors
guides <- endoapply(vectors, FUN = get_guide_seq,
                    upstream_seq = snakemake@params$upstream_seq,
                    downstream_seq = snakemake@params$downstream_seq)

# sort by guide id
guides <- guides[sort(names(guides))]

# get chromosomes for all targeting guides ----------------------------------------------------------

# get indices for different types of guides
guide_ids <- names(guides)
guides_screen <- grepl(guide_ids, pattern = "^chr.+:.+")
guides_nt <- grepl(guide_ids, pattern = "non-targeting_.+")
guides_ctrl <- !Reduce("|", list(guides_screen, guides_nt))

# get chromosomes for screen guides based on id
chrs_screen <- data.frame(guide = guide_ids[guides_screen],
                          chr = sub("^(chr.+):.+", "\\1", guide_ids[guides_screen]),
                          stringsAsFactors = FALSE)

# get target gene names for positive control guides from guide ids
targets_ctrl <- sub("^([[:alnum:]]+).*$", "\\1", guide_ids[guides_ctrl])
chrs_ctrl <- data.frame(guide = guide_ids[guides_ctrl], gene = targets_ctrl,
                        stringsAsFactors = FALSE)
chrs_ctrl[chrs_ctrl$gene == "HS2", "gene"] <- "HBE1"  # HS2 is the enhancer, HBE1 its main target

# get chromosomes for each guide based on target genes and provided target gene annotation
chrs_ctrl <- left_join(chrs_ctrl, gene_chrs, by = "gene")
chrs_ctrl <- select(chrs_ctrl, -gene)

# combine into one data.frame
guide_chrs <- rbind(chrs_ctrl, chrs_screen)

# get binding site(s) for every guide --------------------------------------------------------------

# get BSgenome object
genome <- getBSgenome(snakemake@params$bsgenome)

# get binding sites for all targeting guides
targeting_guides <- structure(guide_ids[!guides_nt], names = guide_ids[!guides_nt])
guide_coords <- bplapply(targeting_guides, FUN = get_guide_coords, guide_seqs = guides,
                         guide_chrs = guide_chrs, genome = genome)
guide_coords <- GRangesList(guide_coords)

# get targets of newly designed screen guides based on id
guide_coords_screen <- guide_coords[grepl(names(guide_coords), pattern = "^chr.+:.+")]
guide_targets_screen <- names(guide_coords_screen) %>% 
  tibble(guide_id = .) %>% 
  mutate(target = sub("^(chr[[:alnum:]]+):([[:digit:]]+)-([[:digit:]]+).+",
                      "\\1-\\2-\\3", guide_id)) %>% 
  separate(target, into = c("chr", "start", "end"), sep = "-", remove = TRUE) %>% 
  makeGRangesListFromDataFrame(split.field = "guide_id")

# only retain binding sites for screen guides that overlap with the respective targets
guide_coords_screen_filt <- endoapply(guide_coords_screen, FUN = filter_guide_coords_target,
                                      guide_targets = guide_targets_screen)

# get target elements for all guides ---------------------------------------------------------------

# convert screen guide coordinates and targets to GRanges object
guide_coords_screen_filt <- unlist(guide_coords_screen_filt)
guide_targets_screen <- unlist(guide_targets_screen)
guide_targets_screen$guide_id <- names(guide_targets_screen)

# infer targets ids for control guides
guide_coords_ctrl <- guide_coords[!grepl(names(guide_coords), pattern = "^chr.+:.+")]
guide_coords_ctrl <- unlist(GRangesList(guide_coords_ctrl))
names(guide_coords_ctrl) <- sub("^([[:alnum:]]+).+", "\\1", names(guide_coords_ctrl))

# infer target regions based on guide coords
guide_coords_ctrl <- split(guide_coords_ctrl, f = names(guide_coords_ctrl))
guide_targets_ctrl <- unlist(endoapply(guide_coords_ctrl, FUN = ctrl_guide_targets))

# combine controls and screen guides and targets
guide_coords <- c(unlist(guide_coords_ctrl), guide_coords_screen_filt)
guide_targets <- c(guide_targets_ctrl, guide_targets_screen)
names(guide_targets) <- NULL
names(guide_coords) <- NULL

# infer target id for every guide
cres <- guide_targets$guide_id
screen_ids <- grepl(cres, pattern = "^chr[[:alnum:]]+:.+")
cres[screen_ids] <- sub("^(chr[[:alnum:]]+:[[:digit:]]+-[[:digit:]]+).+", "\\1", cres[screen_ids])
cres[!screen_ids] <- sub("^([[:alnum:]]+).+", "\\1", cres[!screen_ids])
mcols(guide_targets)[, "target_id"] <- cres

# add target type (enhancer or promoter)
ctrl_enh <- snakemake@params$ctrl_enhancers
type <- if_else(grepl(cres, pattern = "^chr[[:alnum:]]+:.+") | cres %in% ctrl_enh,
                true = "enhancer", false = "promoter")
guide_targets$target_type <- type

# create output data.frame and save to file --------------------------------------------------------

# transform guide coordinates to data.frame
guides_df <- guide_coords %>% 
  as.data.frame() %>% 
  mutate(start = start - 1) %>% 
  select(chr = seqnames, start, end, guide_id = name, strand, spacer = seq)

# transform guide targets to data.frame
guide_targets_df <- guide_targets %>% 
  as.data.frame() %>% 
  mutate(start = start - 1, strand = ".") %>% 
  select(target_chr = seqnames, target_start = start, target_end = end, target_name = target_id,
         target_strand = strand, guide_id)

# add targets to guide data to create main output table
output <- guides_df %>% 
  left_join(guide_targets_df, by = "guide_id") %>% 
  select(chr, start, end, name = guide_id, strand, spacer, target_chr, target_start, target_end,
         target_name, target_strand) %>% 
  arrange(chr, start, end)
  

# add target type column identifying control and enhancer targeting guides
output <- output %>% 
  mutate(target_type = case_when(
    grepl(name, pattern = "^.+-P1.*$") ~ "TSSCtrl",
    grepl(name, pattern = "^.+-[A|B|C|D]$") ~ "enhCtrl",
    TRUE ~ "enh"
  ))

# create bed format guide coordinates
guides_bed <- guides_df %>% 
  mutate(score = 500) %>% 
  select(chr, start, end, guide_id, score, strand, spacer) %>% 
  arrange(chr, start, end)

# create bed format target coordinates
targets_bed <- guide_targets_df %>% 
  mutate(score = 500) %>% 
  select(target_chr, target_start, target_end, target_name, score, target_strand) %>% 
  distinct() %>% 
  arrange(target_chr, target_start, target_end)

# write output to files
write_tsv(output, file = snakemake@output$guide_targets)
write_tsv(guides_bed, file = snakemake@output$guides_bed, col_names = FALSE)
write_tsv(targets_bed, file = snakemake@output$targets_bed, col_names = FALSE)

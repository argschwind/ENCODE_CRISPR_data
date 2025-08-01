## Combine and pre-process HCT116 enhancer-level FLowFISH pipeline results

# save.image("hct116_dl.rda")
# stop()

library(tidyverse)

## Combine data across loci and replicates ---------------------------------------------------------

# input files containing HCT116 data
indir <- "/oak/stanford/groups/engreitz/Projects/Cohesin/220510_TSS_rescaling_rerun/220512_cohesin_fulltiling_tssscaling_ncsat1_nobinall_results/byExperimentRep" 
infiles <- list.files(indir, pattern = "-0-.*FullEnhancerScore.txt", full.names = TRUE)
names(infiles) <- sub(".FullEnhancerScore.txt", "", basename(infiles))

# load all input files
dat <- infiles %>% 
  lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = "id") %>%
  separate(id, into = c("gene", "rep"), sep = "-0-") %>% 
  mutate(rep = as.numeric(sub("Rep1-", "", rep)))

# add perturbation id and name
dat <- dat %>% 
  mutate(PerturbationTargetID = paste0(chr, ":", start, "-", end),
         name = paste0(gene, "|", PerturbationTargetID)) %>% 
  arrange(name, rep)

## Filter replicates -------------------------------------------------------------------------------

# split data by unique tested E-G pairs
dat_split <- split(dat, f = dat$name)

# filter data
filter_data <- function(pair) {
  if (all(pair$Regulated == FALSE)) {
    output <- pair[1, ]
  } else if (all(pair$Regulated == TRUE)) {
    output <- pair[1, ]
  } else {
    output <- tibble()
  }
  return(output)
}

# apply replicate filtering
dat_filt <- lapply(dat_split, FUN = filter_data)
dat_filt <- bind_rows(dat_filt)

## Reformat and save to output files ---------------------------------------------------------------

# compute logFC effect size
dat_filt <- mutate(dat_filt, EffectSize = log(mean / mean.ctrl))

# reformat for output
dat_filt <- dat_filt %>% 
  mutate(strandPerturbationTarget = ".", EffectSize95ConfidenceIntervalLow = NA_real_,
         EffectSize95ConfidenceIntervalHigh = NA_real_, guideSpacerSeq = NA_character_,
         guideSeq = NA_character_, PowerAtEffectSize10 = NA_real_, PowerAtEffectSize25 = NA_real_,
         PowerAtEffectSize50 = NA_real_, ValidConnection = "TRUE", Notes = NA_character_,
         Reference = NA_character_) %>% 
  select(chrom = chr, chromStart = start, chromEnd = end, name, EffectSize,
         strandPerturbationTarget, PerturbationTargetID, EffectSize95ConfidenceIntervalLow,
         EffectSize95ConfidenceIntervalHigh, measuredGeneSymbol = gene, guideSpacerSeq,
         guideSeq, Significant, pValue = p.ttest, pValueAdjusted = fdr.ttest,
         starts_with("PowerAtEffectSize"), ValidConnection, Notes, Reference, Regulated)

# save to output file
write_tsv(dat_filt, file = snakemake@output[[1]])

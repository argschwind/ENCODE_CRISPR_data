## Functions to perform differential gene expression analysis for enhancer TAP-seq

library(dplyr)

#' Censored mean normalization
#' 
#' Apply library size normalization to DGE data based on censored mean, which excludes highly
#' expressed genes.
#' 
#' @param sce A SingleCellExperiment object containing digital gene expression data for all cells
#'   and genes to be normalized.
#' @param percentile Genes of each cell with expression from that percentile upwards will be
#'   excluded to calculate size factors.
#' @param assay Assay to be normalized (default: counts).
normalize_cens_mean <- function(sce, percentile = 0.9, assay = "counts") {
  
  # extract gene expression data
  counts <- assay(sce, assay)
  
  # function to calculate censored umi counts for one cell
  compute_censored_umis <- function(x) {
    sum(x[x <= quantile(x, probs = percentile)]) + 1
  }
  
  # calculate censored total umis per cell
  cens_umis_per_cell <- apply(counts, MARGIN = 2, FUN = compute_censored_umis)
  
  # calculate normalization factors
  norm_factors <- cens_umis_per_cell / mean(cens_umis_per_cell)
  
  # normalize data of each cell based on computed normalization factors
  normcounts <- t(t(counts) / norm_factors)
  
  # add normalized data and normalization factors to sce object
  assay(sce, "normcounts") <- normcounts
  colData(sce)[, "norm_factors"] <- norm_factors
  
  return(sce)
  
}

#' DESeq size factor normalization
#' 
#' Apply library size normalization to DGE data based on size factor calculation methods provided by
#' DESeq2
#' 
#' @param sce A SingleCellExperiment object containing digital gene expression data for all cells
#'   and genes to be normalized.
#' @param percentile Only genes within specified expression quantile will be used to compute size
#'   factors. Default = c(0, 0.9), which excludes to top 10% expressed genes. Set to c(0, 1) to
#'   include all genes for size factor calculation (DESeq2 default).
#' @param norm_genes (optional) Numeric or logical index vector specifying genes to use for size
#'   factor estimation. Use this to compute size factors from e.g. house-keeping genes or known 
#'   unperturbed genes.
#' @param locfun A function to compute a location for a sample (default: median, see
#'   ?estimateSizeFactorsForMatrix for more information)
#' @param type Compute size factors using standard median ("ratio") or positive counts only
#'   ("poscounts"). See ?estimateSizeFactorsForMatrix for more information.
#' @param assay Name of the assay containing counts to be normalized (default: counts).
normalize_deseq <- function(sce, expr_quantile = c(0, 0.9), norm_genes = NULL,
                            locfunc = stats::median, type = c("poscounts", "ratio"),
                            assay = "counts") {
  
  type <- match.arg(type)
  
  # extract counts to normalize
  counts <- assay(sce, assay)
  
  # if no gene list for normalization is provided, get genes within expression quantile range
  if (is.null(norm_genes)) {
    message("Filtering genes based on expression quantile.")
    norm_genes <- get_genes_expr_quantile(counts, expr_quantile = expr_quantile)
  } else {
    message("Using provided genes for size factor estimation.")
  }
  
  # compute size factors
  message("Estimating size factors.")
  size_factors <- DESeq2::estimateSizeFactorsForMatrix(counts, locfunc = locfunc,
                                                       controlGenes = norm_genes, type = type)
  
  # normalize data of each cell based on computed size factors
  normcounts <- t(t(counts) / size_factors)
  
  # add normalized data and normalization factors to sce object
  assay(sce, "normcounts") <- normcounts
  colData(sce)[, "norm_factors"] <- size_factors
  
  return(sce)
  
}

#' Filter Perturb-seq SCE for minimum UMIs per cell
#' 
#' Filter cells in a SCE object for a minimum and maximum number of total UMIs per cell.
#' 
#' @param sce A Perturb-seq SingleCellExperiment object containing UMI counts.
#' @param min_umis,max_umis Minimum and maximum number of UMI counts per cell.
#' @param recompute Should total UMIs per cell be recomputed? Useful if gene filters have been
#'   applied before.
filter_umis_per_cell <- function(sce, min_umis = 0, max_umis = Inf, recompute = FALSE) {
  
  # recompute umis per cell if specified
  if (recompute == TRUE) {
    message("re-comuting total UMIs per cell...")
    colData(sce)[["total_umis"]] <- colSums(assay(sce, "counts"))
  }
  
  # get cells passing umi filter
  total_umis <- colData(sce)[["total_umis"]]
  umi_filter <- total_umis >= min_umis & total_umis <= max_umis
  
  # filter for these cells
  message("Removing ", sum(!umi_filter), " cells based on 'total UMIs per cell' filter.")
  sce <- sce[, umi_filter]
  
}

#' Filter perturbations for minimum cell numbers
#' 
#' Filter perturbations in a Perturb-seq SCE object for a minimum number of cells
#' 
#' @param sce A SingleCellExperiment object containing gene expression data and perturbation
#'   data as alternative experiments (altExp).
#' @param min_cells Minimum number of cells
#' @param pert_level Based on which perturbation level should be filtered? I.e. the name of the
#'   altExp that should be used for filtering
filter_cells_per_pert <- function(sce, min_cells, pert_level) {
  
  # get perturbations with at least 'min_cells' cells
  perts <- altExp(sce, pert_level)
  cells_per_pert <- rowSums(assay(perts, "perts") > 0)
  perts_filter <- cells_per_pert >= min_cells
  
  # filter perturbation altExp based for perturbations passing 'min_cells' filter
  message("Removing ", sum(!perts_filter),
          " perturbations based on 'minimum cells per perturbation' filter.")
  altExp(sce, pert_level) <- perts[perts_filter, ]
  
  return(sce)
  
}

#' Test for differential expression
#' 
#' Test for differential expression across all perturbations. For each perturbation all positive
#' cells are selected and tested against cells carrying non-targeting negative controls.
#' 
#' @param sce A SingleCellExperiment object containing gene expression data and perturbation
#'   data as alternative experiments (altExp)
#' @param pert_level Based on which perturbation level should differential expression test be
#'   performed? I.e. the name of the altExp that should be used as perturbation status matrix.
#' @param max_dist Only consider genes within specified distance from perturbation for differential
#'   expression tests (default: NULL). If NULL all genes are tested against all perturbations.
#' @param de_function Function to used for differential expression testing. Needs to take an SCE
#'   object with at least one column called 'pert' in colData providing perturbation status for a
#'   given perturbation. See code of de_MAST() for an example.
#' @param formula Formula for differential expression model (default: ~pert). Ignored when not
#'   appropriate.
#' @param n_ctrl Specifies how many negative control cells should be randomly drawn
#'   (default: 1000). Set to FALSE to use all non-perturbed cells as negative controls.
#' @param cell_batches (optional) Column name in colData specifying batch for each cell. If
#'   specified control cells are sampled from these batches with equal proportions as perturbed
#'   cells.
#' @param p_adj_method Method to use for multiple testing correction. For more details see
#'   \code{\link[stats]{p.adjust}}.
test_differential_expression <- function(sce, pert_level, max_dist = NULL,
                                         de_function = de_MAST, formula = ~pert, n_ctrl = 5000,
                                         cell_batches = NULL, 
                                         p_adj_method = c("fdr", "holm", "hochberg", "hommel",
                                                          "bonferroni", "BH", "BY", "none")) {
  
  # parse arguments and attach required packages (if needed)
  p_adj_method <- match.arg(p_adj_method)
  
  # check colData
  col_names <- colnames(colData(sce)) 
  if ("pert" %in% col_names) stop("'pert' cannot be a colData name, please rename.", call. = FALSE)
  
  # check that formula contains 'pert'
  if (!"pert" %in% labels(terms(as.formula(formula)))) {
    stop("formula needs to contain 'pert' term to test for perturbation effects.", call. = FALSE)
  }
  
  # check that pert_level is valid
  if (!pert_level %in% altExpNames(sce)) {
    stop("pert_level must be one of the altExp names: ", paste(altExpNames(sce), collapse = ", "),
         call. = FALSE)
  }
  
  # get function to generate input data for one perturbation
  if (is.numeric(n_ctrl)) {
    pert_input_function <- pert_input_sampled
    n_ctrl <- as.integer(n_ctrl)
  } else if (n_ctrl == FALSE) {
    pert_input_function <- pert_input
  } else {
    stop("Invalid 'n_ctrl' argument.", call. = FALSE)
  }
  
  # create vector of perturbations to test
  perts <- structure(rownames(altExp(sce, pert_level)), names = rownames(altExp(sce, pert_level)))
  
  # perform differential gene expression test for each perturbation
  output <- bplapply(perts, FUN = test_de, sce = sce, pert_level = pert_level,
                     cell_batches = cell_batches, pert_input_function = pert_input_function,
                     de_function = de_function, max_dist = max_dist, formula = formula,
                     n_ctrl = n_ctrl)
  
  # convert output into one data.frame
  output <- bind_rows(output, .id = "perturbation")
  
  # correct for multiple testing for all performed tests
  if (method != "LFC") {
    output$pval_adj <- p.adjust(output$pvalue, method = p_adj_method)
    output <- output[order(output$pval_adj), ]
  }
  
  # annotate output with enhancer and gene coordinates and compute distance
  output <- annotate_dist_to_gene(output, sce = sce, pert_level = pert_level)
  
  # add perturbation level used for DE tests
  output$pert_level <- pert_level
  
  # add target type information to output if provided
  pert_data <- as.data.frame(rowData(altExp(sce, pert_level)))
  if ("target_type" %in% colnames(pert_data)) {
    target_type <- pert_data[, c("name", "target_type")]
    output <- left_join(output, target_type, by = c("perturbation" = "name"))
  } else {
    output$target_type <- NA_character_
  }
  
  # add average observed gene expression and number of cells per perturbation
  output <- annotate_cells_and_expr(output, sce = sce, pert_level = pert_level)
  
  return(output)
  
}

# DIFFERENTIAL EXPRESSION TESTS ====================================================================

#' MAST DE function
#' 
#' Perform differential expression analysis between perturbed and control cells using MAST.
#'
#' The first column in colData is expected to provide the perturbation status as a factor with two
#' levels, where the first level is non-perturbed cells. If this is not the case, the function will
#' attempt to convert the first column to a two level factor, with any risks that entails.
#' Any other columns in colData will be used as covariates during modeling, but only the
#' perturbation will be tested for significance.
#'
#' @param pert_object SingleCellExperiment object containing gene expression data and cell groupings
#'   (perturbations) in colData.
#' @param formula Formula to use for differential expression tests. Default: '~ pert', which tests
#'   for an effect of the perturbation status stored in colData as 'pert'.
#' @param pvalue Which test should be used to compute the p-value (default = "hurdle")?
#' @param parallel Should multiple cores be used for fitting (default = FALSE)?
de_MAST <- function(pert_object, formula = ~pert, pvalue = c("hurdle", "cont", "disc"),
                    parallel = FALSE) {
  
  # parse input arguments
  pvalue <- match.arg(pvalue)
  
  # add some row- and colData expected by MAST (not used in DE tests)
  rowData(pert_object) <- cbind(rowData(pert_object), primerid = rownames(pert_object))
  colData(pert_object) <- cbind(colData(pert_object), wellKey = colnames(pert_object))
  
  # coerce pert_object to SingleCellAssay, since MAST requires that as input
  sca <- SceToSingleCellAssay(pert_object)
  
  # fit hurdle model
  zlm_fit <- zlm(as.formula(formula), sca = sca, parallel = parallel)
  
  # calculate log fold changes and confidence intervals using summary function. this will fail if
  # only one gene was tested... in that case return NA
  lfc <- tryCatch(
    withCallingHandlers({
      summary_zlm_fit <- summary(zlm_fit, parallel = parallel)
      summary_dt <- summary_zlm_fit$datatable
      lfc <- summary_dt[contrast == "pert1" & component == "logFC", .(primerid, coef, ci.hi, ci.lo)]
      as.data.frame(lfc)
    }), warning = function(w) {
      message("For perturbation ", pert, ": ", w)
      invokeRestart("muffleWarning")
    }, error = function(e) {
      message("Can't compute logFC: ", e)
      genes <- rownames(pert_object)
      data.frame(primerid = genes, coef = NA_real_, ci.hi = NA_real_, ci.lo = NA_real_)
    })
  
  # perform likelihood ratio test for the perturbation coefficient
  message("Calculating likelihood ratio tests")
  zlm_lr <- lrTest(zlm_fit, "pert")
  pvalues <- zlm_lr[, , "Pr(>Chisq)"]
  if (is.array(pvalues)) {
    pvalues <- data.frame(primerid = names(pvalues[, pvalue]), pvalue = pvalues[, pvalue])
  } else {
    pvalues <- data.frame(primerid = rownames(pert_object), pvalue = pvalues[[pvalue]])
  }
  
  # combine log fold changes and p-values to create output
  output <- merge(lfc, pvalues, by = "primerid")
  colnames(output) <- c("gene", "logFC", "ci_high", "ci_low", "pvalue")
  output <- output[order(output$pvalue), ]
  
  return(output)
  
}


#' DEsingle DE function
#' 
#' The first column in colData is expected to provide the perturbation status as a factor with two
#' levels, where the first level is non-perturbed cells. If this is not the case, the function will
#' attempt to convert the first column to a two level factor, with any risks that entails.
#' Any other columns in colData will ne ignored, because DEsingle currently does not support
#' covariates.
#'
#' @param pert_object SingleCellExperiment object containing gene expression data and cell groupings
#'   in colData. The perturbation to be tested is assumed to be the first column in colData!
#' @param assay Name of the assay used to perform DE tests (default: "normcounts"). Needs to be raw
#'   or size normalized counts (will be rounded to whole numbers).
de_DEsingle <- function(pert_object, assay = "normcounts", ...) {
  
  # extract normalized counts and round to whole numbers
  normcounts <- as.matrix(round(assay(pert_object, assay)))
  
  # get identity for each cell and make sure it's a factor
  groups <- colData(pert_object)[, 1]
  groups <- as.factor(groups)
  levels(groups) <- c(0, 1)
  
  # detect DE genes
  de_results <- DEsingle(counts = normcounts, group = groups)
  
  # reformat output
  output <- data.frame(gene = rownames(de_results), de_results, stringsAsFactors = FALSE,
                       row.names = NULL)
  
  return(output)
  
}

#' LFC DE function
#' 
#' Simple calculation of log fold changes in gene expression between perturbed and control cells.
#' Counts are normalized using censored mean normalization excluding top 10% expressed genes and
#' are log transformed. This is identical to the normalization applied when using MAST for
#' differential expression testing.
#' 
#' @param pert_object SingleCellExperiment object containing gene expression data (counts) and cell
#'   groupings in colData. The perturbation to be tested is assumed to be the first column in
#'   colData!
#' @param assay Name of the assay used to compute log fold change (default: "logcounts").
#' @param pseudocount Pseudocount to be added to transcript counts when calculating average gene
#'   expression.
de_LFC <- function(pert_object, assay = "logcounts", pseudocount = 1, ...) {
  
  # get counts and groups for each cell
  counts <- assay(pert_object, "logcounts")
  groups <- colData(pert_object)[, 1]
  
  # calculate average expression for each gene in both groups
  pert_avg <- rowMeans(counts[, groups == 1] + pseudocount)
  ctrl_avg <- rowMeans(counts[, groups == 0] + pseudocount)
  
  # merge into one data.frame and calculate lfc for each gene
  avg_genex <- data.frame(pert = pert_avg, ctrl = ctrl_avg)
  avg_genex$lfc <- log2(avg_genex$pert) - log2(avg_genex$ctrl)
  
  # reformat output
  output <- data.frame(gene = rownames(avg_genex), avg_genex, stringsAsFactors = FALSE,
                       row.names = NULL)
  
  return(output)
  
}

# HELPER FUNCTIONS =================================================================================

# get genes within expression percentiles based on total expression across all cells
get_genes_expr_quantile <- function(expr, expr_quantile = c(0, 0.9)) {
  
  # get genes within desired expression percentile
  total_expr <- rowSums(expr)
  quants <- quantile(total_expr, probs = expr_quantile)
  gene_filt <- total_expr >= quants[[1]] & total_expr <= quants[[2]]
  
  return(gene_filt)
  
}

# filter a Perturb-seq SCE object for genes within a specified distance of a given perturbation
filt_max_dist_pert <- function(sce, pert_level, pert, max_dist) {
  
  # get genomic coordinates of perturbation pert
  pert_annot <- rowData(altExp(sce, pert_level))
  pert_annot <- makeGRangesFromDataFrame(pert_annot[pert, ])
  
  # only retain data on genes within max_dist from perturbation
  pert_window <- resize(pert_annot, width = max_dist * 2, fix = "center")
  sce <- subsetByOverlaps(sce, pert_window, ignore.strand = TRUE)
  
  return(sce)
  
}

# perform differential gene expression for one perturbation
test_de <- function(pert, sce, pert_level, cell_batches, pert_input_function, max_dist, de_function,
                    formula, n_ctrl) {
  
  # create input object for DE tests
  pert_object <- pert_input_function(pert, sce = sce, pert_level = pert_level,
                                     cell_batches = cell_batches, n_ctrl = n_ctrl)
  
  
  # only retain genes within maximum distance if specified
  if (!is.null(max_dist)) {
    message("Filtering for genes within ", max_dist, " basepairs from perturbation.")
    pert_object <- filt_max_dist_pert(pert_object, pert_level = pert_level, pert = pert,
                                      max_dist = max_dist)
  }
  
  # perform differential gene expression test. warnings and errors get reported and in case of
  # errors NULL is returned for that perturbation.
  output <- tryCatch(
    withCallingHandlers({
      de_function(pert_object, formula = formula)
    }, warning = function(w) {
      message("For perturbation ", pert, ": ", w)
      invokeRestart("muffleWarning")
    }), error = function(e) {
      message("For perturbation ", pert, ": ", e)
      return(NULL)
    })
  
  return(output)
  
}

# generate input for one perturbation without sampling control cells (all non-perturbed cells are
# included as controls)
pert_input <- function(pert, sce, pert_level, ...) {
  
  message("Creating input for perturbation '", pert, "'.")
  
  # get perturbed cells for specified perturbation
  pert_data <- assay(altExp(sce, pert_level), "perts")
  pert_data <- pert_data[pert, ]
  
  # add perturbation status as first column to colData
  pert_status <- DataFrame(pert = as.factor(pert_data[colnames(sce)]))
  colData(sce) <- cbind(pert_status, colData(sce))
  
  # add perturbation id and perturbation level as metadata
  metadata(sce) <- c(metadata(sce), pert_id = pert, pert_level = pert_level)
  
  return(sce)
  
}

# generate input for one perturbation
pert_input_sampled <- function(pert, sce, pert_level, cell_batches, n_ctrl) {
  
  message("Creating input for perturbation '", pert, "' with ", n_ctrl, " sampled control cells.")
  
  # get perturbed cells for specified perturbation
  pert_data <- assay(altExp(sce, pert_level), "perts")
  pert_data <- pert_data[pert, ]
  pert_cells <- pert_data > 0
  
  # randomly draw 'n_ctrl' control cells with same batch distribution as perturbed cells if provided
  if (!is.null(cell_batches)) {
    
    message("Sampling control cells with equal proportions as perturbed cells across '",
            cell_batches, "'.")
    
    # get cell_batches proportions of perturbed cells
    batches <- colData(sce)[cell_batches]
    batches$cell <- rownames(batches)
    pert_cells_prop <- table(batches[pert_cells, cell_batches]) / sum(pert_cells)
    
    # compute the number of control cells per batch to sample with equal proportions
    ctrl_cell_numbers <- round(pert_cells_prop * n_ctrl)
    
    # randomly draw these numbers of control cells from each batch
    non_pert_cells <- batches[!pert_cells, ]
    ctrl_cells <- lapply(names(ctrl_cell_numbers), FUN = function(batch) {
      cells_batch <- non_pert_cells[non_pert_cells[[cell_batches]] == batch, "cell"]
      sample(cells_batch, size = ctrl_cell_numbers[batch], replace = FALSE)
    })
    ctrl_cells <- unlist(ctrl_cells)
    
  } else { 
    
    # simply sample n_ctrl non-perturbed cells if no cell batches are provided
    ctrl_cells <- sample(colnames(sce[, !pert_cells]), size = n_ctrl, replace = FALSE)
    
  }
  
  # extract data for perturbed and sampled control cells
  pert_object <- sce[, c(names(pert_cells)[pert_cells], ctrl_cells)]
  
  # add perturbation status as first column to colData
  pert_status <- DataFrame(pert = as.factor(pert_data[colnames(pert_object)]))
  colData(pert_object) <- cbind(pert_status, colData(pert_object))
  
  # add perturbation id and perturbation level as metadata
  metadata(pert_object) <- c(metadata(pert_object), pert_id = pert, pert_level = pert_level)
  
  return(pert_object)
  
}

# Helper functions to annotate differential expression output --------------------------------------

# annotate differential expression output with number of cells per perturbation and observed average
# gene expression
annotate_cells_and_expr <- function(de_output, sce, pert_level) {
  
  # compute number of cells per perturbation
  perturb_status <- assay(altExp(sce, pert_level), "perts")
  cells_per_pert <- rowSums(perturb_status > 0)
  cells_per_pert <- data.frame(perturbation = names(cells_per_pert), cells = cells_per_pert,
                               stringsAsFactors = FALSE, row.names = NULL)
  
  # compute observed average expression per gene
  avg_expr <- rowMeans(assay(sce, "counts"))
  avg_expr <- data.frame(gene = names(avg_expr), avg_expr = avg_expr, stringsAsFactors = FALSE,
                         row.names = NULL)
  
  # add these to DE output
  de_output <- left_join(de_output, cells_per_pert, by = "perturbation")
  de_output <- left_join(de_output, avg_expr, by = "gene")
  
  return(de_output)
  
}

# annotate differential expression output with perturbation and gene/TSS coordinates as stored in a
# perturb-seq SCE object
annotate_dist_to_gene <- function(de_output, sce, pert_level) {
  
  # get tss coordinates from rowRanges
  gene_coords <- rowRanges(sce)
  gene_coords <- data.frame(gene_chr = as.character(seqnames(gene_coords)),
                            gene_start = start(gene_coords),
                            gene_end = end(gene_coords),
                            gene_strand = as.character(strand(gene_coords)),
                            gene = names(gene_coords),
                            stringsAsFactors = FALSE)
  
  # get perturbation coordinates from alternative experiment rowData
  pert_row_data <- rowData(altExp(sce, pert_level))
  pert_coords <- data.frame(pert_chr = pert_row_data$chr,
                            pert_start = pert_row_data$start,
                            pert_end = pert_row_data$end,
                            perturbation = rownames(pert_row_data),
                            stringsAsFactors = FALSE)
  
  # annotate differential expression output with perturbation, gene/TSS coordinates and distance
  output <- annotate_de_output(de_output, pert_coords = pert_coords, gene_coords = gene_coords)
  
  return(output)
  
}

# annotate differential expression output with perturbation and gene (or TSS) coordinates and
# calculate distance to gene (or TSS)
annotate_de_output <- function(de_output, pert_coords, gene_coords) {
  
  # add perturbation and gene coordinates to de output
  de_output <- left_join(de_output, pert_coords, by = "perturbation")
  de_output <- left_join(de_output, gene_coords, by = "gene")
  
  # compute perturbation center
  de_output$pert_center <- floor((de_output$pert_start + de_output$pert_end) / 2)
  
  # compute distance of perturbation center to gene coordinates
  de_output <- mutate(de_output, distance = case_when(
    pert_chr != gene_chr ~ NA_real_,
    gene_strand == "-" ~ case_when(
      pert_center < gene_start ~ gene_start - pert_center,
      pert_center > gene_end ~ gene_end - pert_center,
      TRUE ~ 0),
    gene_strand == "+" ~ case_when(
      pert_center < gene_start ~ pert_center - gene_start,
      pert_center > gene_end ~ pert_center - gene_end,
      TRUE ~ 0),
    TRUE ~ NA_real_
  ))
  
  return(de_output)
  
}

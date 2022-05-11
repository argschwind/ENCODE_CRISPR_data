## Functions to perform power simulations using specified effect sizes

library(DESeq2)
library(scran)
library(dplyr)
library(Matrix)

#' Fit negative binomial distributions using DESeq2
#' 
#' Fit negative binomial distributions to estimate dispersion for every gene in a
#' SingleCellExperiment object.
fit_negbinom_deseq2 <- function(sce, assay = "counts",
                                size_factors = c("ratio", "poscounts", "iterate", "libsize"),
                                fit_type = c("parametric", "local", "mean"),
                                disp_type = c("dispersion", "dispFit", "dispGeneEst", "dispMAP")) {
  
  size_factors <- match.arg(size_factors)
  fit_type <- match.arg(fit_type)
  disp_type <- match.arg(disp_type)
  
  # check if sce already contains dispersion and raise warning if data will be overwritten
  present_row_data <- colnames(rowData(sce)) %in% c("mean", "dispersion", "disp_outlier_deseq2")
  present_col_data <- colnames(colData(sce)) == "size_factors"
  if (any(present_row_data, present_col_data)) {
    warning("Dispersion data found in sce, will overwrite values", call. = FALSE)
  }
  
  # create DESeq2 object containing count data
  dds <- DESeqDataSetFromMatrix(countData = assay(sce, assay),
                                colData = colData(sce),
                                design = ~ 1)
  
  # compute size factors
  if (size_factors == "libsize") {
    total_umis <- colSums(assay(dds, assay))
    manual_size_factors <- total_umis / mean(total_umis)
    sizeFactors(dds) <- manual_size_factors
  } else {
    dds <- estimateSizeFactors(dds, type = size_factors)
  }
  
  # estimate dispersion
  dds <- estimateDispersions(dds, fitType = fit_type)
  
  # add mean expression, dispersion and cell size factors to rowData and colData of sce
  rowData(sce)[, "mean"] <- rowData(dds)[, "baseMean"]
  rowData(sce)[, "dispersion"] <- rowData(dds)[, disp_type]
  rowData(sce)[, "disp_outlier_deseq2"] <- rowData(dds)[, "dispOutlier"]
  colData(sce)[, "size_factors"] <- sizeFactors(dds)
  
  # store dispersion function in sce
  metadata(sce)[["dispersionFunction"]] <- dispersionFunction(dds)
  
  return(sce)
  
}

#' Differential expression analysis using simulated data sets
#' 
#' Simulate Perturb-seq data sets with a provided effect size and perform differential gene expression
#' analysis as with real data. These results are used to compute the power to detect changes in
#' expression due to CRE perturbations.
#' 
#' @param sce A SingleCellExperiment object containing gene expression data and perturbation
#'   data as alternative experiments (altExp). Needs to contain estimated mean and dispersion for
#'   every gene in rowData.
#' @param effect_size Decrease in gene expression to simulate relative to observed expression (e.g.
#'   0.75 for a 25% decrease in expression).
#' @param pert_level Based on which perturbation level should differential expression test be
#'   performed? I.e. the name of the altExp that should be used as perturbation status matrix.   
#' @param max_dist Only consider genes within specified distance from perturbation for differential
#'   expression tests (default: NULL). If NULL all genes are tested against all perturbations.
#' @param genes_iter Should perturbation effects be simulated one gene at a time (TRUE), or all
#'   genes per perturbation at once (FALSE)? (default: FALSE).
#' @param guide_sd Standard deviation for guide-guide variability simulations. Must be >= 0
#'   (default: 0). 
#' @param center Should average effect size of perturbed cells be centered on specified effect size?
#'   (default: FALSE). Used to make sure that the simulated average effect size does not get shifted
#'   by introduced guide-guide variability.
#' @param rep How many times should data be simulated and differential expression tests be
#'   performed?
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
simulate_diff_expr <- function(sce, effect_size, pert_level, max_dist = NULL, genes_iter = FALSE,
                               guide_sd = 0, center = FALSE, rep = 1,
                               norm = c("real", "sim_nonpert"), de_function = de_MAST,
                               formula = ~ pert, n_ctrl = 5000, cell_batches = NULL) {
  
  # parse input ------------------------------------------------------------------------------------
  
  # parse arguments
  norm <- match.arg(norm)
  
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
  
  # check that guide_sd is valid
  if (!is.numeric(guide_sd) | guide_sd < 0) {
    stop("Invalid 'guide_sd' value. Must be numeric >= 0.", call. = FALSE)
  }
  
  # get required functions -------------------------------------------------------------------------
  
  # get function to generate input data for one perturbation
  if (is.numeric(n_ctrl)) {
    pert_input_function <- pert_input_sampled
    n_ctrl <- as.integer(n_ctrl)
  } else if (n_ctrl == FALSE) {
    pert_input_function <- pert_input
  } else {
    stop("Invalid 'n_ctrl' argument.", call. = FALSE)
  }
  
  # get simulation function based on norm argument
  if (norm == "real") {
    sim_function <- simulate_diff_expr_pert_real
  } else if (norm == "sim_nonpert") {
    sim_function <- simulate_diff_expr_pert_sim
  } else {
    stop("Invalid 'norm' argument.", call. = FALSE)
  }
  
  # perform simulated DE tests ---------------------------------------------------------------------
  
  # create guide_targets data.frame
  guide_targets <- create_guide_targets(sce, pert_level = pert_level)
  
  # create vector of perturbations to test
  perts <- rownames(altExp(sce, pert_level))
  names(perts) <- perts
  
  # simulate Perturb-seq data and perform DE tests for every perturbation
  output <- bplapply(perts, FUN = simulate_diff_expr_rep, sim_function = sim_function, rep = rep,
                     sce = sce, pert_level = pert_level, cell_batches = cell_batches,
                     pert_input_function = pert_input_function, guide_targets = guide_targets, 
                     genes_iter = genes_iter, effect_size = effect_size, guide_sd = guide_sd,
                     center = center, de_function = de_function, max_dist = max_dist,
                     formula = formula, n_ctrl = n_ctrl)
  
  # combine output into one data.frame and rearrange columns
  output <- bind_rows(output, .id = "perturbation")
  output <- relocate(output, iteration, perturbation, gene, perturbed)
  
  return(output)
  
}

#' Repeatedly perform simulated differential expression tests for one perturbation.
#' 
#' Simulate Perturb-seq data and perform differential expression tests a specified number of times.
simulate_diff_expr_rep <- function(pert, sim_function, rep, sce, pert_level, cell_batches,
                                   pert_input_function, guide_targets, genes_iter, effect_size,
                                   guide_sd, center, de_function, max_dist, formula, n_ctrl) {
  
  # repeatedly simulate Perturb-seq data for given perturbation and perform DE tests
  output <- replicate(
    n = rep,
    expr = simulate_diff_expr_try(pert, sim_function = sim_function, sce = sce,
      pert_level = pert_level, cell_batches = cell_batches,
      pert_input_function = pert_input_function, guide_targets = guide_targets,
      genes_iter = genes_iter, effect_size = effect_size, guide_sd = guide_sd, center = center,
      de_function = de_function, max_dist = max_dist, formula = formula, n_ctrl = n_ctrl),
    simplify = FALSE)
  
  # combine output into one data.frame
  output <- bind_rows(output, .id = "iteration")
  output$iteration <- as.integer(output$iteration)
  
  return(output)
  
}

# simulate Perturb-seq data for a given perturbation, perform DE tests and catch any errors or
# warnings. returns 'NULL' if the simulations fail.
simulate_diff_expr_try <- function(pert, sim_function, sce, pert_level, cell_batches,
                                   pert_input_function, guide_targets, genes_iter, effect_size,
                                   guide_sd, center, de_function, max_dist, formula, n_ctrl) {
  
  # try to simulate Perturb-seq data for the given perturbation and catch any errors and warnings
  output <- tryCatch(
    withCallingHandlers({
      sim_function(pert, sce = sce, pert_level = pert_level, cell_batches = cell_batches,
        pert_input_function = pert_input_function, guide_targets = guide_targets,
        genes_iter = genes_iter, effect_size = effect_size, guide_sd = guide_sd, center = center,
        de_function = de_function, max_dist = max_dist, formula = formula, n_ctrl = n_ctrl)
    }, warning = function(w) {
      message("For perturbation ", pert, ": ", w)
      invokeRestart("muffleWarning")
    }), error = function(e){
      message("For perturbation ", pert, ": ", e)
      return(NULL)
    })
  
  return(output)
  
}

## Simulate perturb-seq data and perform DE tests with real data size factors and cell stats -------

# Simulate Perturb-seq data for one perturbation for with an effect of specified effect size on all
# genes within selected maximum distance. Perform differential tests on simulated data using the DE
# function provided by de_function. Perturbation effects can be simulated for all genes at once, or
# one gene at a time, depending on genes_iter parameter.
simulate_diff_expr_pert_real <- function(pert, sce, pert_level, cell_batches, pert_input_function,
                                         guide_targets, genes_iter, effect_size, guide_sd, center,
                                         de_function, max_dist, formula, n_ctrl) {
  
  # create input object for DE tests
  pert_object <- pert_input_function(pert, sce = sce, pert_level = pert_level,
                                     cell_batches = cell_batches, n_ctrl = n_ctrl)
  
  # get perturbation status and gRNA perturbations for all cells
  pert_status <- colData(pert_object)$pert
  grna_perts <- assay(altExp(pert_object, "grna_perts"), "perts")
  
  # get guide ids of guides targeting the perturbed regulatory element
  pert_guides <- guide_targets[guide_targets$target_id == pert, "grna_id"]
  
  # create vector with the gRNA perturbation status for each cell
  grna_pert_status <- create_guide_pert_status(pert_status, grna_perts = grna_perts,
                                               pert_guides = pert_guides)
  
  # only retain genes within maximum distance if specified
  if (!is.null(max_dist)) {
    pert_object <- filt_max_dist_pert(pert_object, pert_level = pert_level, pert = pert,
                                      max_dist = max_dist)
  }
  
  # simulate perturb-seq data and perform differential expression test either all genes at once or
  # one gene at a time
  message("Simulating differential expression.")
  pert_genes <- rownames(pert_object)
  if (genes_iter == FALSE) pert_genes <- list(pert_genes)
  output <- lapply(pert_genes, FUN = simulate_pert_object_real, pert_object = pert_object,
                   effect_size = effect_size, grna_pert_status = grna_pert_status,
                   pert_guides = pert_guides, guide_sd = guide_sd, center = center,
                   pert_status = pert_status, de_function = de_function, formula = formula)
  
  # combine output into one data.frame if simulations were run one perturbed gene at a time
  output <- bind_rows(output)
  
  return(output)
  
}

# simulate expression changes for a set of genes for one perturbation in a pert_object SCE and
# perform differential expression test using real data size factors and cell stats
simulate_pert_object_real <- function(pert_object, pert_genes, effect_size,
                                      grna_pert_status, pert_guides, guide_sd, center,
                                      pert_status, de_function, formula) {
  
  # effect sizes for selected of genes to perturb
  effect_sizes <- structure(rep(1, nrow(pert_object)), names = rownames(pert_object))
  effect_sizes[pert_genes] <- effect_size
  
  # create effect size matrix
  es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                      gene_effect_sizes = effect_sizes, guide_sd = guide_sd)
  
  # set guide variability for non-perturbed guide-gene pairs
  # TODO: implement this in create_effect_size_matrix
  es_mat[, pert_status == 0] <- 1
  es_mat[!rownames(es_mat) %in% pert_genes, ] <- 1
  
  # center effect sizes on specified gene-level effect sizes
  if (center == TRUE) {
    es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status,
                                        gene_effect_sizes = effect_sizes)
  }
  
  # simulate Perturb-seq count data
  sim_object <- sim_tapseq_sce(pert_object, effect_size_mat = es_mat)
  
  # normalize counts using real data normalization factors and log transform
  norm_factors <- colData(sim_object)[, "norm_factors"]
  assay(sim_object, "normcounts") <- t(t(assay(sim_object, "counts")) / norm_factors)
  assay(sim_object, "logcounts") <- log1p(assay(sim_object, "normcounts"))
  
  # perform differential gene expression test
  output <- de_function(sim_object, formula = formula)
  
  # add column labeling genes that were perturbed
  output$perturbed <- output$gene %in% pert_genes  
  
  return(output)
  
}

## Simulate perturb-seq data and perform DE tests with simulated size factors and cell stats -------

# Simulate Perturb-seq data for one perturbation for with an effect of specified effect size on all
# genes within selected maximum distance. Perform differential tests on simulated data using the DE
# function provided by de_function. Perturbation effects can be simulated for all genes at once, or
# one gene at a time, depending on genes_iter parameter.
simulate_diff_expr_pert_sim <- function(pert, sce, pert_level, cell_batches, pert_input_function,
                                        guide_targets, genes_iter, effect_size, guide_sd, center,
                                        de_function, max_dist, formula, n_ctrl) {
  
  # create input object for DE tests
  pert_object <- pert_input_function(pert, sce = sce, pert_level = pert_level,
                                     cell_batches = cell_batches, n_ctrl = n_ctrl)
  
  # get perturbation status and gRNA perturbations for all cells
  pert_status <- colData(pert_object)$pert
  grna_perts <- assay(altExp(pert_object, "grna_perts"), "perts")
  
  # get guide ids of guides targeting the perturbed regulatory element
  pert_guides <- guide_targets[guide_targets$target_id == pert, "grna_id"]
  
  # create vector with the gRNA perturbation status for each cell
  grna_pert_status <- create_guide_pert_status(pert_status, grna_perts = grna_perts,
                                               pert_guides = pert_guides)
  
  # get genes within maximum distance of perturbation to simulate perturbations, else use all genes
  if (!is.null(max_dist)) {
    region_genes <- genes_within_dist(pert_object, pert_level = pert_level, pert = pert,
                                      max_dist = max_dist)
  } else {
    region_genes <- rownames(pert_object)
  }
  
  # simulate perturb-seq data and perform differential expression test either all genes at once or
  # one gene at a time
  message("Simulating differential expression.")
  pert_genes <- region_genes
  if (genes_iter == FALSE) pert_genes <- list(pert_genes)
  output <- lapply(pert_genes, FUN = simulate_pert_object_sim, pert_object = pert_object,
                   region_genes = region_genes, effect_size = effect_size,
                   grna_pert_status = grna_pert_status, pert_guides = pert_guides,
                   guide_sd = guide_sd, center = center, pert_status = pert_status,
                   de_function = de_function, formula = formula)
  
  # combine output into one data.frame if simulations were run one perturbed gene at a time
  output <- bind_rows(output)
  
  return(output)
  
}

# simulate expression changes for a set of genes for one perturbation in a pert_object SCE and
# perform differential expression test using real data size factors and cell stats
simulate_pert_object_sim <- function(pert_object, pert_genes, region_genes, effect_size,
                                     grna_pert_status, pert_guides, guide_sd, center,
                                     pert_status, de_function, formula) {
  
  # effect sizes for selected of genes to perturb
  effect_sizes <- structure(rep(1, nrow(pert_object)), names = rownames(pert_object))
  effect_sizes[pert_genes] <- effect_size
  
  # create effect size matrix
  es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                      gene_effect_sizes = effect_sizes, guide_sd = guide_sd)
  
  # set guide variability for non-perturbed guide-gene pairs
  # TODO: implement this in create_effect_size_matrix
  es_mat[, pert_status == 0] <- 1
  es_mat[!rownames(es_mat) %in% pert_genes, ] <- 1
  
  # center effect sizes on specified gene-level effect sizes
  if (center == TRUE) {
    es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status,
                                        gene_effect_sizes = effect_sizes)
  }
  
  # simulate Perturb-seq count data
  sim_object <- sim_tapseq_sce(pert_object, effect_size_mat = es_mat)
  
  # compute cell stats based on simulated but non-perturbed genes
  ctrl_genes <- !rownames(sim_object) %in% pert_genes
  sim_object$total_umis <- colSums(assay(sim_object[ctrl_genes, ], "counts"))
  sim_object$detected_genes <- colSums(assay(sim_object[ctrl_genes, ], "counts") > 0)
  
  # normalize simulated data using size factors computed on simulated but non-perturbed genes
  sim_object$norm_factors <- normalize_cens_mean(sim_object[ctrl_genes, ])$norm_factors
  assay(sim_object, "normcounts") <- t(t(assay(sim_object, "counts")) / sim_object$norm_factors)
  assay(sim_object, "logcounts") <- log1p(assay(sim_object, "normcounts"))
  
  # perform differential gene expression test
  output <- de_function(sim_object[region_genes, ], formula = formula)
  
  # add column labeling genes that were perturbed
  output$perturbed <- output$gene %in% pert_genes  
  
  return(output)
  
}

## -------------------------------------------------------------------------------------------------

# create a data.frame with gRNAs and their targets
create_guide_targets <- function(sce, pert_level) {
  
  # make gRNAs the targets if pert_level is set to gRNAs, else use targeted elements
  targets <- ifelse(pert_level == "grna_perts", yes = "name", no = "target_name")
  
  # create guide targets data.frame from gRNA annotations
  grnas_annot <- rowData(altExp(sce, "grna_perts"))
  guide_targets <- data.frame(grna_id = grnas_annot[, "name"], target_id = grnas_annot[, targets],
                              stringsAsFactors = FALSE)
  
  return(guide_targets)
  
}

# get genes within a specified distance from a perturbation
genes_within_dist <- function(sce, pert_level, pert, max_dist) {
  
  # get genomic coordinates of perturbation pert
  pert_annot <- rowData(altExp(sce, pert_level))
  pert_annot <- makeGRangesFromDataFrame(pert_annot[pert, ])
  
  # get genes within max_dist from perturbation
  pert_window <- resize(pert_annot, width = max_dist * 2, fix = "center")
  genes <- names(subsetByOverlaps(rowRanges(sce), pert_window, ignore.strand = TRUE))
  
  return(genes)
  
}


## GENERALIZABLE SIMULATION FUNCTIONS ==============================================================

#' Simulate Perturb-seq data based on a SingleCellExperiment object with added mean expression,
#' dispersion and size factors from original real data
#' 
#' @param sce Perturb-seq SingleCellExperiment object
#' @param effect_size_mat Effect size matrix
sim_tapseq_sce <- function(sce, effect_size_mat) {
  
  # simulate Perturb-seq count data with parameters from SCE object
  sim_counts <- simulate_tapseq_counts(gene_means = rowData(sce)[, "mean"],
                                       gene_dispersions = rowData(sce)[, "dispersion"],
                                       cell_size_factors = colData(sce)[, "size_factors"],
                                       effect_size_mat = effect_size_mat)
  
  # convert to SingleCellExperiment object with colData and rowData from sce
  output <- SingleCellExperiment(assays = list(counts = sim_counts), rowData = rowData(sce),
                                 colData = colData(sce))
  
  return(output)
  
}

#' Simulate Perturb-seq data
#' 
#' Simulate Perturb-seq UMI counts for a given number of perturbed cells with a specified effect size on
#' specified genes.
#' 
#' @return A matrix with simulated Perturb-seq UMI counts
simulate_tapseq_counts <- function(gene_means, gene_dispersions, cell_size_factors, effect_size_mat,
                                   gene_ids = names(gene_means),
                                   cell_ids = names(cell_size_factors)) {
  
  # number of genes and cells
  n_genes <- length(gene_means)
  n_cells <- length(cell_size_factors)
  
  # make mu matrix for simulation
  mu <- matrix(rep(gene_means, n_cells), ncol = n_cells)
  mu <- sweep(mu, 2, cell_size_factors, "*")  # add cell-to-cell variability based on size factors
  
  # inject perturbation effects by element-wise product of mu and effect_size_mat
  mu <- mu * effect_size_mat
  
  # simulate counts
  Matrix(rnbinom(n_cells * n_genes, mu = mu, size = 1 / gene_dispersions), ncol = n_cells,
         dimnames = list(gene_ids, cell_ids))
  
}

# guide-guide variability specific functions -------------------------------------------------------

# function to randomly pick 1 expressed guide per cell
sample_guide <- function(pert_status) {
  
  # randomly pick one expressed guide
  apply(pert_status, MARGIN = 2, FUN = function(x) {
    out <- structure(rlang::rep_along(names(x), x = 0), names = names(x))
    guides <- which(x > 0)
    sampled_guide <- guides[sample(length(guides), size = 1)]
    out[sampled_guide] <- 1
    return(out)
  })
  
}

# create vector with the gRNA perturbation status for each cell
create_guide_pert_status <- function(pert_status, grna_perts, pert_guides) {
  
  # get gRNA perturbations for perturbed cells
  grnas_pert_cells <- grna_perts[pert_guides, pert_status == 1]
  
  # if >1 guides target the given perturbation convert guide perturbation matrix into vector with
  # unique guide-level perturbation status for each cell. If a cell expresses multiple guides, one
  # is randomly selected
  if (!is.null(nrow(grnas_pert_cells))) {
    grnas_pert_cells <- convert_pert_mat_to_vector(grnas_pert_cells)
  }
  
  # get gRNA perturbations for control cells and also convert these into a vector if needed
  grnas_ctrl_cells <- grna_perts[!rownames(grna_perts) %in% pert_guides, pert_status == 0]
  if (!is.null(nrow(grnas_ctrl_cells))) {
    grnas_ctrl_cells <- convert_pert_mat_to_vector(grnas_ctrl_cells)
  }
  
  # adjust control gRNA status so that they 'come after' targeting gRNAs
  ctrl_perts <- grnas_ctrl_cells > 0
  grnas_ctrl_cells[ctrl_perts] <- grnas_ctrl_cells[ctrl_perts] + max(grnas_pert_cells)
  
  # combine perturbed and control perturbation status vectors
  c(grnas_pert_cells, grnas_ctrl_cells)
  
}

# convert a perturbation status matrix to a vector with a unique status for every perturbation. If a
# cell has >1 perturbations, one is chosen randomly
convert_pert_mat_to_vector <- function(pert_mat) {
  
  # randomly pick one perturbation if a cell has > 1
  pert_mat_sampled <- sample_guide(pert_mat)
  
  # remove perturbation that do not occur anymore
  pert_mat_sampled <- pert_mat_sampled[rowSums(pert_mat_sampled) > 0, ]
  
  # create unique perturbation status for each perturbation and transform to vector
  pert_mat_unique_status <- sweep(pert_mat_sampled, 1, seq_len(nrow(pert_mat_sampled)), "*")
  colSums(pert_mat_unique_status)
  
}

# function to create an effect size matrix with added guide-guide variability in effect size
create_effect_size_matrix <- function(grna_pert_status, pert_guides, gene_effect_sizes, guide_sd) {
  
  # randomly draw effect size variation of guides on every gene
  n_pert_guides <- length(pert_guides)
  n_ctrl_guides <- max(grna_pert_status) - n_pert_guides
  guide_effect_sizes_pert <- vapply(gene_effect_sizes, FUN = rnorm, n = n_pert_guides,
                                    sd = guide_sd, FUN.VALUE = numeric(n_pert_guides))
  guide_effect_sizes_ctrl <- vapply(rlang::rep_along(gene_effect_sizes, 1), FUN = rnorm,
                                    n = n_ctrl_guides, sd = guide_sd,
                                    FUN.VALUE = numeric(n_ctrl_guides))
  guide_effect_sizes <- rbind(guide_effect_sizes_pert, guide_effect_sizes_ctrl)
  
  # set negative guide effect sizes to 0
  guide_effect_sizes[guide_effect_sizes < 0] <- 0
  
  # add row with no effect for non-perturbed cells
  guide_effect_sizes <- rbind(1, guide_effect_sizes)
  
  # pick correct effect sizes for every cell based on it's gRNA perturbation status
  es_mat <- t(guide_effect_sizes[grna_pert_status + 1, ])
  colnames(es_mat) <- names(grna_pert_status)
  
  return(es_mat)
  
}

# function to center effect size matrix with guide variability so that the average effect size per
# gene corresponds to a specified effect size
center_effect_size_matrix <- function(effect_size_mat, pert_status, gene_effect_sizes) {
  
  # get mean effect size for every gene for perturbed and control cells
  mean_es_pert <- rowMeans(effect_size_mat[, pert_status == 1])
  mean_es_ctrl <- rowMeans(effect_size_mat[, pert_status == 0])
  
  # compute required shift to center guide effect sizes on the specified effect sizes
  pert_shift <- gene_effect_sizes - mean_es_pert
  ctrl_shift <- 1 - mean_es_ctrl
  
  # center guide-level effect sizes on specified effect sizes
  effect_size_mat[, pert_status == 1] <- effect_size_mat[, pert_status == 1] + pert_shift
  effect_size_mat[, pert_status == 0] <- effect_size_mat[, pert_status == 0] + ctrl_shift
  
  # set negative guide effect sizes due to shift to 0
  effect_size_mat[effect_size_mat < 0] <- 0
  
  return(effect_size_mat)
  
}

## DEPRECATED FUNCTIONS ============================================================================

# create effect size matrix without guide variability
# create_effect_size_matrix <- function(pert_status, gene_effect_sizes) {
#   
#   # convert pert_status to logical vector
#   pert_status <- pert_status == 1
#   
#   # create matrix with specified effect sizes for genes in perturbed cells
#   es_mat <- matrix(1, ncol = length(pert_status), nrow = length(gene_effect_sizes),
#                    dimnames = list(names(gene_effect_sizes), names(pert_status)))
#   es_mat[, pert_status] <- sweep(es_mat[, pert_status], 1, gene_effect_sizes, "*")
#   
#   return(es_mat)
#   
# }

#' Fit negative binomial distributions using monocle
#' 
#' Fit negative binomial distributions to estimate dispersion for every gene in a
#' SingleCellExperiment object.
# fit_negbinom_monocle <- function(sce, disp_type = c("dispersion_empirical", "dispersion_fit"),
#                                  remove_outliers = TRUE) {
#   
#   disp_type <- match.arg(disp_type)
#   
#   # check if sce already contains dispersion and raise warning if data will be overwritten
#   present_row_data <- colnames(rowData(sce)) %in% c("mean", "dispersion")
#   present_col_data <- colnames(colData(sce)) == "size_factors"
#   if (any(present_row_data, present_col_data)) {
#     warning("Dispersion data found in sce, will overwrite values", call. = FALSE)
#   }
#   
#   # create monocle object containing transcript counts
#   genes <- data.frame(gene_short_name = rownames(sce), row.names = rownames(sce),
#                       stringsAsFactors = FALSE)
#   cds <- newCellDataSet(cellData = as.matrix(assay(sce, "counts")),
#                         featureData = new("AnnotatedDataFrame", data = genes))
#   
#   # estimate size factors and dispersion
#   cds <- estimateSizeFactors(cds)
#   cds <- estimateDispersions(cds, remove_outliers = remove_outliers)
#   
#   # add mean expression, dispersion and cell size factors to rowData and colData of sce
#   rowData(sce)[, "mean"] <- dispersionTable(cds)[, "mean_expression"]
#   rowData(sce)[, "dispersion"] <- dispersionTable(cds)[, disp_type]
#   colData(sce)[, "size_factors"] <- pData(cds)[, "Size_Factor"]
#   
#   # store dispersion function in sce
#   metadata(sce)[["dispersionFunction"]] <- cds@dispFitInfo$blind$disp_func
#   
#   return(sce)
#   
# }

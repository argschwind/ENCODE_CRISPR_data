
#' Test for trans-acting differential expression effects
test_trans_diff_expression <- function(sce, perts, pert_level, sample_genes = NULL,
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
  
  # perform differential gene expression test for each perturbation
  output <- bplapply(perts, FUN = test_de_trans, sce = sce, pert_level = pert_level,
                     cell_batches = cell_batches, pert_input_function = pert_input_function,
                     de_function = de_function, formula = formula, sample_genes = sample_genes,
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

# perform differential gene expression for trans-effects
test_de_trans <- function(pert, sce, pert_level, cell_batches, pert_input_function, de_function,
                          formula, n_ctrl, sample_genes) {
  
  # create input object for DE tests
  pert_object <- pert_input_function(pert, sce = sce, pert_level = pert_level,
                                     cell_batches = cell_batches, n_ctrl = n_ctrl)
  
  # only retain genes on other chromosomes than the perturbation
  pert_object <- filter_trans_genes(pert_object, pert = pert, pert_level = pert_level,
                                    sample_genes = sample_genes)
  
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

# filter sce object for trans genes for a give perturbation
filter_trans_genes <- function(sce, pert, pert_level, sample_genes) {
  
  # get chromosome of perturbation
  pert_chr <- rowData(altExp(sce, pert_level))[pert, "chr"]
  
  # get all trans genes
  genes <- rowRanges(sce)
  trans_genes <- names(genes[seqnames(genes) != pert_chr])
  
  # sample to a specific number of genes if specified
  if (!is.null(sample_genes)) {
    trans_genes <- sample(trans_genes, size = sample_genes, replace = FALSE)
  }
  
  # filter sce for trans genes only
  sce <- sce[trans_genes, ]
  
  return(sce)
  
}

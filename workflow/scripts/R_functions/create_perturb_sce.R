## Functions to create Perturb-seq SCE objects

library(SingleCellExperiment)
library(Matrix)
library(tidyr)
library(dplyr)

#' Create a perturb-seq SCE object
#' 
#' Create an SingleCellObject containing expression data, gene annotations and perturbation status
#' for each cell.
#' 
#' @param expr Matrix containing transcript counts for each cell (rows = genes, columns = cells).
#' @param annot Named GRanges object with genome annotations for genes (rows) in expr. It's strongly
#'   suggested that annot is named the same way as the rows in expr (e.g. gene ids) to ensure
#'   correct association of annotations with genes in expression data.
#' @param pert_status A data.frame (or tibble) containing gRNAs detected in each cell (long format).
#' @param grna_annot A data.frame (or tibble) with gRNA coordinates (columns 1-7) and optional
#'   target element annotations (columns 8-13). Used to extract genomic coordinates for
#'   perturbations. If target columns are provided, also used to create target-level perturbation
#'   status, e.g. if a cell carries one or more gRNAs targeting specified cis-regulatory elements.
#' @param targets_pert_name How should the target-level perturbations be called?
#'   (default: "target_perts")
#' @value An SCE object.
create_pert_sce <- function(expr, annot, pert_status, grna_annot,
                            targets_pert_name = "target_perts") {
  
  # check grna_annot format
  validate_grna_annot(grna_annot)
  
  # check for cells and gRNAs consistency
  if (length(setdiff(pert_status$cell, colnames(expr))) > 0) {
    stop("Not all cells in pert_status found in expr.", call. = FALSE)
  }
  
  if (length(setdiff(pert_status$grna, grna_annot$name)) > 0) {
    stop("Not all gRNAs in pert_status are found in targets.", call. = FALSE)
  }
  
  # process feature annotations --------------------------------------------------------------------
  
  # if annot does not have names raise warning and abort if annot contains more features than genes.
  # if annot has names, check that names of annot uniquely identify genes
  if (is.null(names(annot))) {
    if (length(annot) == nrow(expr)) {
      warning("annot does not have names. Assuming annot is in order with rows in expr.",
              call. = FALSE)
    } else {
      stop("annot does not have names and differs in length from rows in expr. ", 
           "Cannot assign annotations to genes in expr.", call. = FALSE)
    }
  } else if (all(rownames(expr) %in% names(annot)) == FALSE) {
    stop("Not all expr rownames found in annot names. Cannot assign annotations to genes in expr.",
         call. = FALSE)
  } else if (any(table(names(annot)[names(annot) %in% rownames(expr)]) > 1)) {
    stop("Duplicated names in annot for genes in expr. Cannot assign annotations to genes in expr.",
         call. = FALSE)
  } else {
    annot <- annot[rownames(expr)]
  }
  
  # compute number of UMIs and detected genes per cell ---------------------------------------------
  
  # compute total umis and number of detected genes per cell
  umis_per_cell <- colSums(expr)
  genes_per_cell <- colSums(expr > 0)
  
  # create DataFrame with cell stats
  cell_stats <- DataFrame(total_umis = umis_per_cell,
                          detected_genes = genes_per_cell,
                          row.names = names(umis_per_cell) )
  
  # create perturbation status SummarizedExperiment objects ----------------------------------------
  
  # create gRNA-level perturbation matrix
  grna_pert_status <- create_pert_matrix(pert_status, expr = expr, pert_ids = "grna")
  
  # create perturbation status SE object for gRNA-level perturbations
  grna_perts_se <- create_pert_se_obj(grna_pert_status, pert_annot = grna_annot)
  
  # create target-level SE object, if gRNA targets are provided
  target_cols <- c("target_chr", "target_start", "target_end", "target_name", "target_strand")
  if (all(target_cols %in% colnames(grna_annot)) == TRUE) {
    
    # add gRNA targets for each gRNA to pert_status
    pert_status <- left_join(pert_status, grna_annot[, c("name", "target_name")],
                             by = c("grna" = "name"))
    
    # create target-level perturbation matrix
    target_pert_status <- create_pert_matrix(pert_status, expr = expr, pert_ids = "target_name")
    
    # get unique targets and rename mandatory columns
    target_annot <- distinct(grna_annot[, seq(from = 7, to = ncol(grna_annot), by = 1)])
    cols_idx <- colnames(target_annot) %in% target_cols
    colnames(target_annot)[cols_idx] <- sub("target_", "", colnames(target_annot)[cols_idx])
    
    # create perturbation status SE object for target-level perturbations
    target_perts_se <- create_pert_se_obj(target_pert_status, pert_annot = target_annot)
    
  } else {
    target_perts_se <- NULL
  }
  
  # create SingleCellExperiment output -------------------------------------------------------------
  
  # create SingleCellExperiment object from expression, column data and gene annotations
  sce <- SingleCellExperiment(assays = list(counts = expr), colData = cell_stats, rowRanges = annot)
  
  # add perturbations to gene expression SCE object as alternative experiment(s)
  altExp(sce, e = "grna_perts") <- grna_perts_se
  altExp(sce, e = targets_pert_name) <- target_perts_se
  
  return(sce)
  
}


# HELPER FUNCTIONS =================================================================================

# create perturbation status matrix from a pert_status table and a gene expression matrix
create_pert_matrix <- function(pert_status, expr, pert_ids = c("grna", "target_name"),
                               pert_counts = FALSE) {
  
  # parse pert_ids argument
  pert_ids <- match.arg(pert_ids)
  
  # select cells and perturbations and remove NAs if there are any
  pert_status <- pert_status[, c("cell", pert_ids)]
  if (any(is.na(pert_status[[pert_ids]]))) {
    warning("NAs found in ", pert_ids, ". These perturbations will be dropped!", call. = FALSE)
    pert_status <- drop_na(pert_status)
  }
  
  # remove duplicated perturbations if perturbations shouldn't be counted per cell
  if (pert_counts == FALSE) pert_status <- distinct(pert_status)
  
  # get unique cells and target ids
  cells <- unique(colnames(expr))
  target_ids <- unique(pert_status[[pert_ids]])
  
  # create row and column indices identifying the coordinates of each cell - target combination in 
  # the sparse matrix to be created
  pert_status$row_index <- match(pert_status[[pert_ids]], target_ids)
  pert_status$col_index <- match(pert_status$cell, cells)
  
  # create sparse perturbation status matrix
  pert_status_matrix <- sparseMatrix(
    i = pert_status$row_index, 
    j = pert_status$col_index,
    x = 1L, 
    dims = c(length(target_ids), length(cells)),
    dimnames = list(target_ids, cells)
  )
  
  return(pert_status_matrix)
  
}

# create a perturbation status SummarizedExperiment object
create_pert_se_obj <- function(pert_status_mat, pert_annot) {
  
  # create DataFrames from pert_annot to use as rowData in SummarizedExperiment objects
  pert_annot <- DataFrame(pert_annot, row.names = pert_annot$name)
  
  # create SummarizedExperiment objects containing perturbation status and annotations
  pert_se_obj <- SummarizedExperiment(assays = list(perts = pert_status_mat),
                                      rowData = pert_annot[rownames(pert_status_mat), ])
  
  return(pert_se_obj)
  
}

# validate grna_annot input
validate_grna_annot <- function(grna_annot) {
  
  # allowed column name sets
  names_6cols <- c("chr", "start", "end", "name", "strand", "spacer")
  names_7cols <- c(names_6cols, "target_type")
  names_11cols <- c(names_6cols, paste0("target_", names_6cols[1:5]))
  names_12cols <- c(names_11cols, "target_type")
  
  # create list with all column name options
  col_options <- list(names_6cols, names_7cols, names_11cols, names_12cols)
  
  # check if grna column names are valid
  cols <- colnames(grna_annot)
  cols_valid <- vapply(col_options, FUN = identical, cols, FUN.VALUE = logical(1))
  
  # raise error if column names do not match one of the options
  if (all(cols_valid == FALSE)) {
    stop("Invalid column names/order in grna_annot: [", paste(cols, collapse = ", "), "]",
         call. = FALSE)
  }
  
}

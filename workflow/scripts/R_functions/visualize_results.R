## Functions to visualize perturbation - target pairs

library(ggplot2)
library(dplyr)
library(cowplot)
library(SingleCellExperiment)

# plot one or more perturbation-gene pairs from differential expression results
plot_pairs <- function(pairs, sce, assay = "counts", n_ctrl = NULL, plot = TRUE) {
  
  # raise warning if number of pairs is large
  if (nrow(pairs) > 50) {
    warning("Large number of pairs (", nrow(pairs), ")")
  }
  
  # get column for gene, perturbation, perturbation level and adjusted p-value (since apply
  # converts rows to vectors). this prevents errors if column order in pairs changes
  gene_col = which(colnames(pairs) == "gene")
  pert_col = which(colnames(pairs) == "perturbation")
  pert_level_col = which(colnames(pairs) == "pert_level")
  pval_col = which(colnames(pairs) == "pval_adj")
  
  # plot all pairs
  p <- apply(pairs, MARGIN = 1, FUN = plot_pair_apply, sce = sce, assay = assay, n_ctrl = n_ctrl,
             gene_col = gene_col, pert_col = pert_col, pert_level_col = pert_level_col,
             pval_col = pval_col)
  
  # arrange plots (or return list of plots)
  if (plot == TRUE) {
    plot_grid(plotlist = p)
  } else {
    return(p)
  }

}

# basic function to plot gene expression of a gene (y) as function perturbation status (x). this
# function returns a ggplot object, which allows adding or overwriting elements of the plot
plot_pert_gene_pair <- function(gene, pert, sce, pert_level, assay = "counts", n_ctrl = NULL,
                                title = paste(gene, pert, sep = " | ")) {
  
  # get perturbation status and gene expression
  pert_status <- assay(altExp(sce, pert_level), assay = "perts")[pert, ]
  gene_expr <- assay(sce, assay)[gene, ]
  
  # get perturbed and control cells
  pert_cells <- names(pert_status[pert_status > 0])
  ctrl_cells <- names(pert_status[pert_status == 0])
  
  # sample n_ctrl control cells (if specified)
  if (!is.null(n_ctrl)) ctrl_cells <- sample(ctrl_cells, size = n_ctrl)

  # create vector with all cells to plot
  cells <- c(pert_cells, ctrl_cells)
  
  # extract data for selected cells and create data frame
  dat <- tibble(
    cell = cells,
    pert_status = factor(pert_status[cells], levels = sort(unique(pert_status))),
    gene_expr = gene_expr[cells]
  )
  
  # create colors for perturbation status
  pert_colors <- c("gray55", rep("firebrick2", n_distinct(dat$pert_status) - 1))

  # plot gene expression as a function of perturbation status
  ggplot(dat, aes(x = pert_status, y = gene_expr, group = pert_status, color = pert_status)) +
    geom_jitter(alpha = 0.5) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +
    labs(title = title, x = "Perturbation status", y = paste0("Gene expression (", assay, ")")) +
    scale_color_manual(values = pert_colors) +
    theme_bw() +
    theme(legend.position = "none")
  
}

## HELPER FUNCTIONS ================================================================================

# plot one perturbation - gene pair from differential expression results. used to apply plotting to
# many pairs
plot_pair_apply <- function(pair, sce, assay, n_ctrl, gene_col, pert_col, pert_level_col, pval_col) {
  
  # get gene, perturbation ids and adjusted p-value from pair
  gene <- pair[[gene_col]]
  pert <- pair[[pert_col]]
  pert_level <- pair[[pert_level_col]]
  pval <- pair[[pval_col]]
  
  # create title
  title = paste0(gene, " | ", pert, " (p.adj = ", pval, ")")
  
  # create plot
  plot_pert_gene_pair(gene, pert = pert, sce = sce, pert_level = pert_level, assay = assay,
                      n_ctrl = n_ctrl, title = title)
  
}

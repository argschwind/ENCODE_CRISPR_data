
# required packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)
library(cowplot)

# estimate guide variability using a linear model and fit normal distribution to output
estimate_guide_variability <- function(per_target, per_guide, guide_targets, cells_per_guide = 25,
                                       effect_size = c("logFC", "log2FC", "pctChange"),
                                       pct_change_range = c(-1, 0),
                                       return_plots = FALSE) {
  
  # only retain per-guide results with a specified minimum of cells
  per_guide <- filter(per_guide, cells >= cells_per_guide)
  
  # merge per-target and per-guide differential expression results
  merged <- combine_de_results(per_target = per_target, per_guide = per_guide,
                               guide_targets = guide_targets, effect_size = effect_size)
  
  # filter pairs for per-target percentage change range 
  merged <- merged %>% 
    filter(target_effect_size >= pct_change_range[[1]], target_effect_size <= pct_change_range[[2]])
  
  # fit linear model to estimate guide variability
  guide_var <- fit_guide_var_model(merged)
  
  # fit normal distribution and extract fitted parameters
  distr_fit <- fitdistr(guide_var$guide_var, densfun = "normal")
  distr_param <- distr_fit$estimate
  
  # plot per-target vs per-guide variability and fitted regression model
  p1 <- plot_guide_var(guide_var)
  
  # plot gRNA variability distribution
  p2 <- plot_guide_var_distr(guide_var, distr_param = distr_param)
  
  # create output
  output <- list(guide_variability = guide_var,
                 distribution = as_tibble_row(distr_param))
  
  # add plots to output if specified, else plot them
  if (return_plots == TRUE) {
    output$plots <- list(guide_variability = p1, distribution = p2)
  } else {
    plot_grid(p1, p2, rel_widths = c(4, 3))
  }
  
  return(output)
         
}

# combine per-target and per-guide differential expression results
combine_de_results <- function(per_target, per_guide, guide_targets,
                               effect_size = c("logFC", "log2FC", "pctChange")) {
  
  effect_size <- match.arg(effect_size)
  
  # convert effect sizes to percent change (if needed)
  if (effect_size != "pctChange") {
    base <- switch(effect_size, "logFC" = exp(1), "log2FC" = 2)
    per_target$effect_size <- logFC_to_pctChange(per_target$effect_size, base = base)
    per_guide$effect_size  <- logFC_to_pctChange(per_guide$effect_size, base = base)
  }
  
  # filter out any cases with NA effect size
  per_target <- filter(per_target, !is.na(effect_size))
  per_guide <- filter(per_guide, !is.na(effect_size))
  
  # select required column to merge de results
  per_target_merge <- dplyr::select(per_target, target = perturbation, gene,
                                    target_effect_size = effect_size, pval_adj)
  
  per_guide_merge <- dplyr::select(per_guide, guide = perturbation, gene,
                                   guide_effect_size = effect_size)
  
  # add target to per guide results
  per_guide_merge <- left_join(per_guide_merge,
                               dplyr::select(guide_targets, guide = name, target = target_name),
                               by = "guide")
  
  # merge per target and per guide results and rearrange columns
  merged <- inner_join(per_target_merge, per_guide_merge, by = c("target", "gene"))
  merged <- relocate(merged, guide, target, gene, pval_adj, target_effect_size, guide_effect_size)
  
  return(merged)
  
}

# fit linear model to model per-guide effect size as a function of per-target effect size. regress 
# out the per-target effect from per-guide effect size to estimate variability due to different
# guides
fit_guide_var_model <- function(merged) {
  
  # fit regression model
  model <- lm(guide_effect_size ~ target_effect_size, data = merged)
  
  # extract residuals (guide effect size normalized for per-target effect) and add to merged data
  guide_var <- residuals(model)
  merged$guide_var <- guide_var
  
  return(merged)
  
}

# convert log-fold change to percentage change
logFC_to_pctChange <- function(logFC, base = exp(1)) {
  pct_change <- -(1 - base ^ logFC)
  return(pct_change)
}

# plot guide variabilities
plot_guide_var <- function(guide_var) {
  ggplot(guide_var, aes(x = target_effect_size, y = guide_effect_size)) +
    geom_point(alpha = 0.25) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    labs(title = "Per-target vs per-guide effect size",
         x = "per-target effect size (pct change)",
         y = "per-guide effect size (pct change)") +
    theme_bw() +
    theme(text = element_text(size = 10))
}

# plot guide variability distribution
plot_guide_var_distr <- function(guide_var, distr_param) {
  ggplot(guide_var, aes(guide_var)) +
    geom_histogram(aes(y =..density..), bins = 100) +
    labs(title = "Guide variability distribution",
         x = "Guide effect size normalized for target effect") +
    stat_function(fun = dnorm, args = list(mean = distr_param[["mean"]], sd = distr_param[["sd"]]),
                  color = "firebrick3", lwd = 1) +
    theme_bw() +
    theme(text = element_text(size = 10))
}

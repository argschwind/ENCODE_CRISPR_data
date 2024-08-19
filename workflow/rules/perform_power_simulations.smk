## Rules to perform empirical power simulations for Perturb-seq experiment

# this workflow allows splitting of datasets by chromosome for easy parallelization. whether a
# sample should be split by chromosomes is specified in the config file via the 'split_by_chr'
# entry. rule ordering is used to apply chromosome splitting for samples that are set for splitting.
# wildcard constraints are used to define which samples should be split.

# set rule order so that rules splitting a sample by chromosome are preferred over other rules
ruleorder: compute_power_chrs > compute_power

# get samples that should be split by chromsome and create wildcards constrain string
split_samples = list(dict(filter(lambda x: x[1] == True, config["split_by_chr"].items())).keys())
split_samples_wildcards = "|".join(split_samples)

# Rules for power simulations ----------------------------------------------------------------------

# fit negative binomial distribution to estimate dispersions
rule fit_negbinom_distr:
  input: "resources/{sample}/perturb_sce.rds"
  output: temp("resources/{sample}/perturb_sce_disp.rds")
  params:
    umis_per_cell = lambda wildcards: config["diff_expr"]["umis_per_cell"][wildcards.sample],
    remove_genes = lambda wildcards: config["power_simulations"]["remove_genes"][wildcards.sample],
    size_factors = config["power_simulations"]["size_factors"],
    fit_type = config["power_simulations"]["fit_type"]
  log: "results/{sample}/logs/power_sim/fit_negbinom_distr.log"
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "64G",
    runtime = "6h"
  script:
    "../scripts/fit_negbinom_distr.R"

# perform power simulations by submitting one iteration at a time
rule perform_power_simulations:
  input: "resources/{sample}/perturb_sce_disp.rds"
  output:
    temp("results/{sample}/power_sim/rep{rep}_output_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz")
  params:
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = lambda wildcards: config["diff_expr"]["max_dist"][wildcards.sample],
    formula = config["diff_expr"]["formula"],
    n_ctrl = lambda wildcards: config["diff_expr"]["n_ctrl"][wildcards.sample],
    norm = lambda wildcards: config["power_simulations"]["norm"][wildcards.sample],
    cell_batches = lambda wildcards: config["diff_expr"]["cell_batches"][wildcards.sample],
    genes_iter = False
  log: "results/{sample}/logs/power_sim/power_sim_rep{rep}_{effect}_{sd}gStd_{method}_{strategy}.log.gz"
  threads: config["power_simulations"]["threads"]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "48G",
    runtime = "24h"
  script:
   "../scripts/power_simulations.R"

# compute power
rule compute_power:
  input:
    expand("results/{{sample}}/power_sim/rep{rep}_output_{{effect}}_{{sd}}gStd_{{method}}_{{strategy}}.tsv.gz",
      rep = range(1, config["power_simulations"]["rep"] + 1))
  output:
    "results/{sample}/power_sim/power_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz"
  params:
    p_adj_method = config["diff_expr"]["padj_method"],
    pval_threshold = config["diff_expr"]["padj_threshold"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/compute_power.R"

# Rules for simulations split by chromosome --------------------------------------------------------

# fit negative binomial distribution to estimate dispersions
rule fit_negbinom_distr_chr:
  input: "resources/{sample}/perturb_sce.{chr}.rds"
  output: temp("resources/{sample}/perturb_sce_disp.{chr}.rds")
  params:
    umis_per_cell = lambda wildcards: config["diff_expr"]["umis_per_cell"][wildcards.sample],
    remove_genes = lambda wildcards: config["power_simulations"]["remove_genes"][wildcards.sample],
    size_factors = config["power_simulations"]["size_factors"],
    fit_type = config["power_simulations"]["fit_type"]
  log: "results/{sample}/logs/power_sim/fit_negbinom_distr_{chr}.log"
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "32G",
    runtime = "5h"
  script:
    "../scripts/fit_negbinom_distr.R"

# perform power simulations by submitting one iteration at a time
rule perform_power_simulations_chr:
  input: "resources/{sample}/perturb_sce_disp.{chr}.rds"
  output:
    temp("results/{sample}/power_sim/{chr}/{chr}_rep{rep}_output_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz")
  params:
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = lambda wildcards: config["diff_expr"]["max_dist"][wildcards.sample],
    formula = config["diff_expr"]["formula"],
    n_ctrl = lambda wildcards: config["diff_expr"]["n_ctrl"][wildcards.sample],
    norm = lambda wildcards: config["power_simulations"]["norm"][wildcards.sample],
    cell_batches = lambda wildcards: config["diff_expr"]["cell_batches"][wildcards.sample],
    genes_iter = False
  log: "results/{sample}/logs/power_sim/{chr}/power_sim_{chr}_rep{rep}_{effect}_{sd}gStd_{method}_{strategy}.log.gz"
  threads: config["power_simulations"]["threads"]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "96G",
    runtime = "24h"
  script:
   "../scripts/power_simulations.R"

# compute power
rule compute_power_chrs:
  input:
    expand("results/{{sample}}/power_sim/{chr}/{chr}_rep{rep}_output_{{effect}}_{{sd}}gStd_{{method}}_{{strategy}}.tsv.gz",
      chr = ["chr" + str(i)  for i in [*range(1, 23), "X"]],
      rep = range(1, config["power_simulations"]["rep"] + 1))
  output:
    "results/{sample}/power_sim/power_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz"
  wildcard_constraints:
    sample = split_samples_wildcards
  params:
    p_adj_method = config["diff_expr"]["padj_method"],
    pval_threshold = config["diff_expr"]["padj_threshold"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/compute_power.R"
    
# Additional rules ---------------------------------------------------------------------------------

# estimate guide-guide variability from differential expression results
rule estimate_guide_variability:
  input:
    per_target = "results/{sample}/diff_expr/output_{method}_perCRE.tsv.gz",
    per_guide = "results/{sample}/diff_expr/output_{method}_perGRNA.tsv.gz",
    guide_targets = "resources/{sample}/guide_targets.tsv"
  output: 
    guide_var = "results/{sample}/guide_var/guide_variability_{method}.tsv",
    distr_fit = "results/{sample}/guide_var/guide_variability_distribution_{method}.tsv",
    plots = "results/{sample}/guide_var/guide_variability_{method}.pdf"
  params:
    effect_size_col = "logFC",
    fdr_sig = 0.05,
    cells_per_guide = 25,
    pct_change_range = [-0.5, -0.1]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/estimate_guide_variability.R"

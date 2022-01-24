## Rules to perform empirical power simulations for Perturb-seq experiment

# extimate guide-guide variability from differential expression results
rule estimate_guide_variability:
  input:
    per_target = "results/{sample}/diff_expr/output_{method}_perCRE.csv.gz",
    per_guide = "results/{sample}/diff_expr/output_{method}_perGRNA.csv.gz",
    guide_targets = "resources/{sample}/guide_targets.tsv"
  output: 
    guide_var = "results/{sample}/guide_variability_{method}.csv",
    distr_fit = "results/{sample}/guide_variability_distribution_{method}.csv",
    plots = "results/{sample}/guide_variability_{method}.pdf"
  params:
    effect_size_col = "logFC",
    fdr_sig = 0.05,
    fdr_nonsig = 0.1,
    cells_per_guide = 25
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/estimate_guide_variability.R"

# fit negative binomial distribution to estimate dispersions
rule fit_negbinom_distr:
  input: "resources/{sample}/perturb_sce.{chr}.rds"
  output: temp("resources/{sample}/perturb_sce_disp.{chr}.rds")
  params:
    umis_per_cell = config["diff_expr"]["umis_per_cell"]["Gasperini"],
    remove_genes = config["power_simulations"]["remove_genes"],
    size_factors = "poscounts",
    fit_type = config["power_simulations"]["fit_type"]
  log: "results/{sample}/logs/power_sim/fit_negbinom_distr_{chr}.log"
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/fit_negbinom_distr.R"

# perform power simulations by submitting one iteration at a time
rule perform_power_simulations:
  input: "resources/{sample}/perturb_sce_disp.{chr}.rds"
  output:
    temp("results/{sample}/power_sim/{chr}_rep{rep}_output_{effect}_{sd}gStd_{method}_{strategy}.csv.gz")
  params:
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    cell_batches = None,
    pert_genes = config["power_simulations"]["pert_genes"]
  log: "results/{sample}/logs/power_sim/power_sim_{chr}_rep{rep}_{effect}_{sd}gStd_{method}_{strategy}.log.gz"
  threads: config["power_simulations"]["threads"]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "24G",
    time = "24:00:00"
  script:
   "../scripts/power_simulations.R"
   
# compute power
rule compute_power:
  input:
    expand("results/{{sample}}/power_sim/{chr}_rep{rep}_output_{{effect}}_{{sd}}gStd_{{method}}_{{strategy}}.csv.gz",
      chr = ["chr" + str(i)  for i in [*range(1, 23), "X"]],
      rep = range(1, config["power_simulations"]["rep"] + 1))
  output:
    "results/{sample}/power_sim/power_{effect}_{sd}gStd_{method}_{strategy}.csv.gz"
  params:
    fdr = 0.05,
    seed = 210928
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/compute_power.R"

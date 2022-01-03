## Rules to perform empirical power simulations for Perturb-seq experiment

# perform power simulations for one chromosome
rule perform_power_simulations:
  input: "resources/{sample}/perturb_sce.{chr}.rds"
  output:
    temp("results/{sample}/power_sim/{chr}_output_{effect}_{sd}gStd_{method}_{strategy}.csv.gz")
  params:
    umis_per_cell = config["diff_expr"]["umis_per_cell"]["Gasperini"],
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    cell_batches = None,
    remove_genes = config["power_simulations"]["remove_genes"],
    size_factors = "poscounts",
    fit_type = config["power_simulations"]["fit_type"],
    rep = config["power_simulations"]["rep"],
    pert_genes = config["power_simulations"]["pert_genes"],
    seed = 20211217
  log: "results/{sample}/logs/power_sim/power_sim_{chr}_{effect}_{sd}gStd_{method}_{strategy}.log.gz"
  threads: config["power_simulations"]["threads"]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "128G",
    time = "72:00:00"
  script:
   "../scripts/power_simulations.R"


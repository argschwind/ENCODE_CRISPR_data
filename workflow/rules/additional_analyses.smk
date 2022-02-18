## Rules for additional analyses concerning simulation strategies

# simulations using different levels of guide variability ------------------------------------------

rule guide_var_simulations_realCovars:
  input: "resources/{sample}/perturb_sce_disp.{chr}.rds"
  output:
    temp("results/{sample}/additional_analyses/guide_var_simulations/realCovars/{chr}_rep{rep}_output_{effect}_{sd}gStd_{method}_{strategy}.csv.gz")
  params:
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    cell_batches = None,
    genes_iter = False,
    norm = "real"
  log: "results/{sample}/guide_var_simulations/logs/power_sim_{chr}_rep{rep}_{effect}_{sd}gStd_{method}_{strategy}.log.gz"
  threads: config["power_simulations"]["threads"]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "4G",
    time = "00:30:00"
  script:
   "../scripts/power_simulations.R"
   
rule guide_var_simulations_simCovars:
  input: "resources/{sample}/perturb_sce_disp.{chr}.rds"
  output:
    temp("results/{sample}/additional_analyses/guide_var_simulations/simCovars/{chr}_rep{rep}_output_{effect}_{sd}gStd_{method}_{strategy}.csv.gz")
  params:
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    cell_batches = None,
    genes_iter = False,
    norm = "sim_nonpert"
  log: "results/{sample}/guide_var_simulations/logs/power_sim_{chr}_rep{rep}_{effect}_{sd}gStd_{method}_{strategy}.log.gz"
  threads: config["power_simulations"]["threads"]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "4G",
    time = "00:30:00"
  script:
   "../scripts/power_simulations.R"

rule guide_var_analysis:
  input:
    expand("results/Gasperini2019/additional_analyses/guide_var_simulations/{covars}/{chr}_rep{rep}_output_{effect}_{sd}gStd_MAST_perCRE.csv.gz",
           covars = ["realCovars", "simCovars"], chr = ["chr21"], effect = [0.25],
           rep = range(1, 51), sd = [0, 0.06, 0.1, 0.2, 0.5])
  output: "results/Gasperini2019/additional_analyses/guide_var_simulations/guide_var_analysis.html"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/gasperini_dataset/guide_var_analysis.Rmd"

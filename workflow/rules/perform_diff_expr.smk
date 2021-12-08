## Rules to perform differential expression on a Perturb-seq experiment to e.g. map CREs to target 
## genes

# extract data for one chromosome
rule extract_chromosome:
  input: "resources/{sample}/perturb_sce.rds"
  output: temp("resources/{sample}/perturb_sce.{chr}.rds")
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "8G"
  script:
    "../scripts/extract_chrom_from_sce.R"

# perform differential expression tests for one chromosome
rule perform_de_tests:
  input: "resources/{sample}/perturb_sce.{chr}.rds"
  output:
    temp("results/{sample}/diff_expr/{chr}_output_{method}_{strategy}.csv.gz")
  params:
    umis_per_cell = config["diff_expr"]["umis_per_cell"]["Gasperini"],
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    cell_batches = None,
    seed = 20210928
  log: "results/{sample}/logs/diff_expr_{chr}_{method}_{strategy}.log"
  threads: config["diff_expr"]["threads"]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "38G",
    time = "3:00:00"
  script:
    "../scripts/differential_expression.R"
    
# run DE tests for all chromosomes and combine into one file
rule combine_de_results:
  input:
    expand("results/{{sample}}/diff_expr/{chr}_output_{{method}}_{{strategy}}.csv.gz",
      chr = ["chr" + str(i) for i in [*range(1, 23), "X"]])
  output:
    "results/{sample}/diff_expr/output_{method}_{strategy}.csv.gz"
  params:
    padj_method = config["diff_expr"]["padj_method"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/combine_diff_expr_tests.R"

## Rules to combine and format differential expression and power analysis output

# compile output files in ENCODE format
rule create_encode_dataset:
  input:
    diff_expr = "results/{sample}/diff_expr/output_{method}_{strategy}.csv.gz",
    power_sim = expand("results/{{sample}}/power_sim/power_{effect}_{{sd}}gStd_{{method}}_{{strategy}}.csv.gz",
                  effect = config["power_simulations"]["effect_sizes"]),
    annot = "resources/gencode.v26lift37.annotation.gtf.gz",
    guide_targets = "resources/{sample}/guide_targets.tsv"
  output: "results/{sample}/full_data_{sd}gStd_{method}_{strategy}.tsv.gz"
  params:
    gene_ids = "ensembl",
    tss_ctrl_tag = "PosCtrl",
    tss_min_dist = 1000,
    padj_threshold = config["diff_expr"]["padj_threshold"],
    cell_type = lambda wildcards: config["metadata"][wildcards.sample]["cell_type"],
    reference = lambda wildcards: config["metadata"][wildcards.sample]["reference"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_output_dataset.R"

# plot power analysis results
rule power_analysis:
  input:
    "results/{sample}/full_data_{sd}gStd_{method}_{strategy}.tsv.gz"
  output:
    "results/{sample}/power_analysis_{sd}gStd_{method}_{strategy}.html"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/power_analysis.Rmd"

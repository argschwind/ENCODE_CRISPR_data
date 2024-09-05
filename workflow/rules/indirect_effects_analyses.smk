## Estimate indirect (trans-acting) effects in CRISPR data and filter based on expected proportion
## of indirect effects

# split whole all perturbations into 10 batches to run in parallel
rule create_pert_batches:
  input: "resources/{sample}/perturb_sce.rds"
  output: temp("results/{sample}/trans_effects/pert_batches_{strategy}.tsv")
  params:
    n_batches = lambda wildcards: config["trans_effects"]["n_batches"][wildcards.sample]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "32G"
  script:
    "../scripts/indirect_effects/create_perturbation_batches.R"

# perform differential expression tests for trans-acting effects for one batch of perturbations
rule perform_trans_effect_de_tests_batch:
  input: 
    sce = "resources/{sample}/perturb_sce.rds",
    pert_batches = "results/{sample}/trans_effects/pert_batches_{strategy}.tsv"
  output: temp("results/{sample}/trans_effects/batch{batch}_output_{method}_{strategy}.tsv.gz")
  params:
    umis_per_cell = lambda wildcards: config["diff_expr"]["umis_per_cell"][wildcards.sample],
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    sample_genes = lambda wildcards: config["trans_effects"]["sample_genes"][wildcards.sample],
    formula = config["diff_expr"]["formula"],
    n_ctrl = lambda wildcards: config["diff_expr"]["n_ctrl"][wildcards.sample],
    cell_batches = lambda wildcards: config["diff_expr"]["cell_batches"][wildcards.sample],
    p_adj_method = config["diff_expr"]["padj_method"],
    seed = 20210928
  log: "results/{sample}/logs/trans_effects/batch{batch}_trans_diff_expr_{method}_{strategy}.log"
  threads: config["diff_expr"]["threads"]
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "84G",
    runtime = "24h"
  script:
    "../scripts/indirect_effects/perform_trans_effects_de_tests.R"

# perform all trans-effect test for Gasperini2019 dataset
rule trans_effects_gasperini:
  input:
    cis_results = "results/Gasperini2019/diff_expr/output_MAST_perCRE.tsv.gz",
    trans_results = expand("results/Gasperini2019/trans_effects/batch{batch}_output_MAST_perCRE.tsv.gz",
                           batch = [1, 2, 3, 4, 5, 6, 8, 9, 10])
  output: "results/Gasperini2019/trans_effects/output_trans_effects_MAST_perCRE.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/indirect_effects/combine_trans_effects_de_results.R"
    
# plot expected proportion of indirect effects against proportion of effects in cis analysis
rule plot_indirect_effects:
  input:
    trans_results = "results/{sample}/trans_effects/output_trans_effects_MAST_perCRE.tsv.gz",
    cis_results = "results/{sample}/diff_expr/output_MAST_perCRE.tsv.gz",
    encode_results = "results/ENCODE/ENCODE_{sample}_{sd}gStd_MAST_perCRE_{genome}.tsv.gz"
  output: "results/{sample}/trans_effects/indirect_effects_{sd}gStd_MAST_perCRE_{genome}.{dist}kb_filter.html"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/indirect_effects/plot_indirect_effects.Rmd"  
  
# filter CRISPR dataset for potential indirect effects
rule filter_indirect_effects_gasperini:
  input: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_Gasperini2019_0.13gStd_0.8pwrAt15effect_GRCh38.tsv.gz"
  output: 
    filt = "results/Gasperini2019/trans_effects/EPCrisprBenchmark_Gasperini2019_0.13gStd_0.8pwrAt15effect_GRCh38.{dist}kb_hardFilter.tsv.gz",
    flip = "results/Gasperini2019/trans_effects/EPCrisprBenchmark_Gasperini2019_0.13gStd_0.8pwrAt15effect_GRCh38.{dist}kb_flipFilter.tsv.gz",
    prop_filt = "results/Gasperini2019/trans_effects/EPCrisprBenchmark_Gasperini2019_0.13gStd_0.8pwrAt15effect_GRCh38.{dist}kb_propFilter.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/indirect_effects/filter_indirect_effects.R"  

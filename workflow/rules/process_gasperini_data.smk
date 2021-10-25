
# download Gasperini et al. 2019 data
rule download_gasperini_cds:
  output: temp("resources/Gasperini2019/GSE120861_at_scale_screen.cds.rds")
  conda: "../envs/r_process_crispr_data.yml"
  params:
    url = config["download_urls"]["gasperini_cds"]
  shell:
    "wget -O {output}.gz {params.url} ; gunzip {output}.gz"
    
rule download_gasperini_guide_targets:
  output: "resources/Gasperini2019/GSE120861_gene_gRNAgroup_pair_table.at_scale.txt.gz"
  conda: "../envs/r_process_crispr_data.yml"
  params:
    url = config["download_urls"]["gasperini_guide_tagets"]
  shell:
    "wget -O {output} {params.url}"
    
# download annotations
rule download_gencode_annotations:
  output: "resources/gencode.v26lift37.annotation.gtf.gz"
  conda: "../envs/r_process_crispr_data.yml"
  params:
    url = config["download_urls"]["gencode_v26lift37"]
  shell:
    "wget -O {output} {params.url}"
    
# convert to summarized experiment object like for TAP-seq data, with added gene annotation
rule create_sce_gasperini:
  input:
    sce = "resources/Gasperini2019/GSE120861_at_scale_screen.cds.rds",
    annot = "resources/gencode.v26lift37.annotation.gtf.gz",
    guide_targets = "resources/Gasperini2019/GSE120861_gene_gRNAgroup_pair_table.at_scale.txt.gz"
  output:
    "resources/Gasperini2019/GSE120861_at_scale_screen.sce.rds"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_sce_gasperini.R"

# extract data for one chromosome
rule extract_chromosome:
  input: "resources/Gasperini2019/GSE120861_at_scale_screen.sce.rds"
  output: temp("resources/Gasperini2019/GSE120861_at_scale_screen.sce.{chr}.rds")
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "8G"
  script:
    "../scripts/extract_chrom_from_sce.R"

# perform differential expression tests for one chromosome
rule perform_de_tests:
  input: "resources/Gasperini2019/GSE120861_at_scale_screen.sce.{chr}.rds"
  output:
    temp("results/Gasperini2019/diff_expr/{chr}_output_{method}_{strategy}.csv.gz")
  params:
    umis_per_cell = config["diff_expr"]["umis_per_cell"]["Gasperini"],
    min_cells = lambda wildcards: config["diff_expr"]["min_cells"][wildcards.strategy],
    max_dist = config["diff_expr"]["max_dist"],
    formula = config["diff_expr"]["formula"],
    n_ctrl = config["diff_expr"]["n_ctrl"],
    cell_batches = None,
    seed = 20210928
  log: "results/Gasperini2019/logs/diff_expr_{chr}_{method}_{strategy}.log"
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
    expand("results/Gasperini2019/diff_expr/{chr}_output_{{method}}_{{strategy}}.csv.gz",
      chr = ["chr" + str(i) for i in [*range(1, 23), "X"]])
  output:
    "results/Gasperini2019/diff_expr/output_{method}_{strategy}.csv.gz"
  params:
    padj_method = config["diff_expr"]["padj_method"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/combine_diff_expr_tests.R"

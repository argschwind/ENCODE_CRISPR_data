## Rules to create worflow input and output

# download UCSC hg19 to hg38 liftover chain file
rule download_chain_file:
  output: "resources/hg19ToHg38.over.chain.gz"
  params:
    url = config["download_urls"]["liftover_chain"]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# download gencode annotations
rule download_gencode_annotations:
  output: "resources/{annot}.annotation.gtf.gz"
  params:
    url = lambda wildcards: config["download_urls"][wildcards.annot]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"
    
# create output data by combining differential expression and power simulation results
rule create_output:
  input:
    diff_expr = "results/{sample}/diff_expr/output_{method}_{strategy}.tsv.gz",
    power_sim = expand("results/{{sample}}/power_sim/power_{effect}_{{sd}}gStd_{{method}}_{{strategy}}.tsv.gz",
                  effect = config["power_simulations"]["effect_sizes"])
  output: "results/{sample}/output_{sd}gStd_{method}_{strategy}.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_output.R"
    
# plot differential expression results
rule diff_expr_results:
  input:
    "results/{sample}/diff_expr/output_{method}_{strategy}.tsv.gz"
  output:
    "results/{sample}/diff_expr_{method}_{strategy}.html"
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "8G"
  script:
    "../scripts/diff_expr_results.Rmd"
    
# plot power analysis results
rule power_analysis:
  input:
    "results/{sample}/output_{sd}gStd_{method}_{strategy}.tsv.gz"
  output:
    "results/{sample}/power_analysis_{sd}gStd_{method}_{strategy}.html"
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "8G"  
  script:
    "../scripts/power_analysis.Rmd"

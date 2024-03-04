
# process HCT116 results and reformat to ENCODE CRISPR format
rule process_hct116_results:
  input:
    results = "resources/HCT116/hct116_data/240227_cohesin_noAux_combined_analysis.csv",
    tss = "resources/HCT116/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed",
    annot = "resources/gencode.v29.annotation.gtf.gz"
  output: "results/ENCODE/ENCODE_HCT116_FlowFISH_GRCh38.tsv.gz"
  params:
    tss_min_dist = config["encode_datasets"]["dist_to_TSS"][0],
    padj_threshold = config["diff_expr"]["padj_threshold"],
    reference = "Guckelberger et al., inPrep"
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/hct116_dataset/process_hct116_results.R"
    
# filter EP benchmarking datasets for distance to TSS and minimum power
rule filter_hct116_results:
  input: "results/ENCODE/ENCODE_HCT116_FlowFISH_GRCh38.tsv.gz"
  output: temp("results/ENCODE/EPCrisprBenchmark/ENCODE_HCT116_FlowFISH_{pwr}pwrAt{es}effect_GRCh38.tsv.gz")
  params:
    tss_to_dist = config["encode_datasets"]["dist_to_TSS"],
    remove_filtered_pairs = False
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/filter_crispr_dataset.R"

# convert ENCODE format files to EPBenchmarking format files
rule create_hct116_ep_benchmarking_dataset:
  input: "results/ENCODE/EPCrisprBenchmark/ENCODE_HCT116_FlowFISH_{pwr}pwrAt{es}effect_GRCh38.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_HCT116_FlowFISH_{pwr}pwrAt{es}effect_GRCh38.tsv.gz"
  params:
    effect_size = "pctChange",
    min_pct_change = None,
    cell_type = "HCT116"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"  

## Combine enhancer-level FlowFISH pipeline output and reformat to EP benchmarking format ----------

rule combine_hct116_replicates:
  output: "results/HCT116/hct116_hg19.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/hct116_dataset/combine_hct116_replicates.R"
  
rule liftover_hct116_enhancers:
  input:
    results = "results/HCT116/hct116_hg19.tsv.gz",
    chain = "resources/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/HCT116/enh_hct116_hg19.bed",
    hg38 = "results/HCT116/enh_hct116_hg38.bed",
    unlifted = "results/HCT116/enh_hct116_unlifted.bed"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "zcat {input.results} | "
    """awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $1":"$2"-"$3, 0, "."}}' | """
    "sort | uniq > {output.hg19} ; "
    "liftOver {output.hg19} {input.chain} {output.hg38} {output.unlifted}"
    
rule liftover_hct116_dataset:
  input:
    results = "results/HCT116/hct116_hg19.tsv.gz",
    enh = "results/HCT116/enh_hct116_hg38.bed",
    tss = "resources/HCT116/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed",
    annot = "resources/gencode.v29.annotation.gtf.gz",
  output: "results/HCT116/EPCrisprBenchmark_HCT116_GRCh38.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/hct116_dataset/liftover_hct116_data.R"

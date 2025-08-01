## Reformat HCT116 CRISPRi FlowFISH results to ENCODE and EP benchmarking CRISPR formats

# process HCT116 results and reformat to ENCODE CRISPR format
rule process_hct116_results:
  input:
    pos = "resources/HCT116/hct116_data/240227_cohesin_noAux_combined_analysis_significant.csv",
    neg = "resources/HCT116/hct116_data/240227_cohesin_noAux_combined_analysis_negatives.csv", 
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

# convert ENCODE format files to EPBenchmarking format files
rule create_hct116_ep_benchmarking_dataset:
  input: "results/ENCODE/ENCODE_HCT116_FlowFISH_GRCh38.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_HCT116_FlowFISH_GRCh38.tsv.gz"
  params:
    effect_size = "pctChange",
    min_pct_change = None,
    cell_type = "HCT116"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"

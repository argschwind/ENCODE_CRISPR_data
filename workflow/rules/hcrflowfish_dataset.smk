## Rules to download and reformat ENCODE HCR-FlowFISH data

# process metadata file containing download URLs
rule process_metadata:
  input: "resources/HCRFlowFish/hcrflowfish_metadata_240123.tsv"
  output: "resources/HCRFlowFish/hcrflowfish_files.tsv"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/hcrflowfish_dataset/parse_metadata.R"

# download all HCR FlowFish files and create ENCODE format output
rule download_processed_hcr_flowfish:
  input:
    urls = "resources/HCRFlowFish/hcrflowfish_files.tsv",
    encode_format = "resources/HCRFlowFish/encode_cre-gene_crispr_format_may2021.tsv",
    cres = "resources/DNase/ENCFF205FNC/ENCFF205FNC_candidate_cres.bed",
    annot = "resources/gencode.v29.annotation.gtf.gz",
    training_crispr_data = "/oak/stanford/groups/engreitz/Projects/Benchmarking/CRISPR_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz"
  output: temp("results/ENCODE/ENCODE_HCRFlowFISH_perCRE.tsv.gz")
  params:
    reference = "Reilly et al., 2021",
    ignore_txs = ["ENST00000380252.6", "ENST00000292896.3", "ENST00000380237.5"],
    tss_min_dist = config["encode_datasets"]["dist_to_TSS"][0],
    padj_threshold = config["diff_expr"]["padj_threshold"]
  threads: 5
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/hcrflowfish_dataset/download_files.R"
    
# filter for cis E-G pairs only
rule filter_hcr_flowfish:
  input: "results/ENCODE/ENCODE_HCRFlowFISH_perCRE.tsv.gz"
  output: temp("results/ENCODE/ENCODE_HCRFlowFISH_perCRE_filt.tsv.gz")
  params:
    tss_to_dist = [1000, 1e6]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/hcrflowfish_dataset/filter_hcr_flowfish.R"

# create EP benchmarking format file
rule create_ep_benchmarking_file:
  input: "results/ENCODE/ENCODE_HCRFlowFISH_perCRE_filt.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_HCRFlowFISH_perCRE_GRCh38.tsv.gz"
  params:
     effect_size = "logFC",
     cell_type = "K562"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"

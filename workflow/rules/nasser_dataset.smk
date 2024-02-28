## Rules to download and format the K562 CRISPR dataset from Nasser et al., 2021 for EP benchmarking

## Download and reformat data into ENCODE CIRSPR format --------------------------------------------

# download Nasser et al., 2021 data
rule download_nasser:
  output: temp("results/ENCODE/Nasser2021_crisprData_{set}.tsv")
  params:
    url = lambda wildcards: config["download_urls"]["Nasser2021"][wildcards.set]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"
    
# download ABC RefSeq TSS annotations
rule download_refseq_tss:
  output: "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed"
  params:
    url = config["download_urls"]["RefSeq_tss"]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# reformat into ENCODE format
rule reformat_nasser:
  input: 
    data = "results/ENCODE/Nasser2021_crisprData_{set}.tsv",
    tss = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed"
  output: temp("results/ENCODE/ENCODE_Nasser2021_{set}_hg19.tsv.gz")
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/nasser_dataset/format_nasser_data.R"

## Lift elements and TSSs over to hg38 -------------------------------------------------------------
    
# lift enhancer coordinates from hg19 to hg38 using UCSC's liftOver software    
rule liftover_nasser_enhancers:
  input:
    results = "results/ENCODE/ENCODE_Nasser2021_{set}_hg19.tsv.gz",
    chain = "resources/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/ENCODE/liftover/Nasser2021/enh_{set}_hg19.bed",
    hg38 = "results/ENCODE/liftover/Nasser2021/enh_{set}_hg38.bed",
    unlifted = "results/ENCODE/liftover/Nasser2021/enh_{set}_unlifted.bed"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "zcat {input.results} | "
    """awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $1":"$2"-"$3, 0, "."}}' | """
    "sort | uniq > {output.hg19} ; "
    "liftOver {output.hg19} {input.chain} {output.hg38} {output.unlifted}"
    
# liftover Nasser2021 data from hg19 to hg38
rule liftover_nasser_dataset:
  input:
    results = "results/ENCODE/ENCODE_Nasser2021_{set}_hg19.tsv.gz",
    enh_hg38 = "results/ENCODE/liftover/Nasser2021/enh_{set}_hg38.bed",
    annot_hg38 = "resources/gencode.v29.annotation.gtf.gz"
  output: temp("results/ENCODE/ENCODE_Nasser2021_{set}_GRCh38.tsv.gz")
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/encode_datasets/liftover_crispr_dataset.R"
    
## Create EPBenchmarking CRISPR data files ---------------------------------------------------------

# convert K562 data in ENCODE format to EPBenchmarking format file
rule create_nasser_ep_benchmarking_dataset:
  input: "results/ENCODE/ENCODE_Nasser2021_K562_GRCh38.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_Nasser2021_K562_GRCh38.tsv.gz"
  params:
     effect_size = "pctChange",
     min_pct_change = None,
     cell_type = "K562"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"
    
# convert data on all cell types to EPBenchmarking format file
rule create_nasser_allCellTypes_ep_benchmarking_dataset:
  input: "results/ENCODE/ENCODE_Nasser2021_allCellTypes_GRCh38.tsv.gz"
  output: temp("results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_Nasser2021_allCellTypes_GRCh38.tsv.gz")
  params:
     effect_size = "pctChange",
     min_pct_change = None,
     cell_type_col = "Notes"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"
    
# subset Nasser et al 2021 dataset into different cell types and convert to EPBenchmarking format
rule subset_nasser_ep_benchmarking_dataset:
  input: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_Nasser2021_allCellTypes_GRCh38.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_Nasser2021_{set}_GRCh38.tsv.gz"
  params:
    subset = lambda wildcards: config["Nasser2021"]["subsets"][wildcards.set]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/nasser_dataset/split_nasser_ep_benchmarking_dataset.R"

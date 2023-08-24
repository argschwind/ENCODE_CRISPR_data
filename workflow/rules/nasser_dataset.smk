## Rules to download and format the K562 CRISPR dataset from Nasser et al., 2021 for EP benchmarking

## Download and reformat data into ENCODE CIRSPR format --------------------------------------------

# download Nasser et al., 2021 data
rule download_nasser:
  output: temp("results/ENCODE/Nasser2021_data.tsv")
  params:
    url = config["download_urls"]["Nasser2021"]
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
    data = "results/ENCODE/Nasser2021_data.tsv",
    tss = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS500bp.bed"
  output: "results/ENCODE/ENCODE_Nasser2021_hg19.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/nasser_dataset/format_nasser_data.R"

## Lift elements and TSSs over to hg38 -------------------------------------------------------------
    
# lift enhancer coordinates from hg19 to hg38 using UCSC's liftOver software    
rule liftover_nasser_enhancers:
  input:
    results = "results/ENCODE/ENCODE_Nasser2021_hg19.tsv.gz",
    chain = "resources/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/ENCODE/liftover/Nasser2021/enh_hg19.bed",
    hg38 = "results/ENCODE/liftover/Nasser2021/enh_hg38.bed",
    unlifted = "results/ENCODE/liftover/Nasser2021/enh_unlifted.bed"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "zcat {input.results} | "
    """awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $1":"$2"-"$3, 0, "."}}' | """
    "sort | uniq > {output.hg19} ; "
    "liftOver {output.hg19} {input.chain} {output.hg38} {output.unlifted}"
    
# liftover Nasser2021 data from hg19 to hg38
rule liftover_nasser_dataset:
  input:
    results = "results/ENCODE/ENCODE_Nasser2021_hg19.tsv.gz",
    enh_hg38 = "results/ENCODE/liftover/Nasser2021/enh_hg38.bed",
    annot_hg38 = "resources/gencode.v29.annotation.gtf.gz"
  output: "results/ENCODE/ENCODE_Nasser2021_GRCh38.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/liftover_crispr_dataset.R"
    
## Create EPBenchmarking CRISPR data file ----------------------------------------------------------

# convert ENCODE format files to EPBenchmarking format files
rule create_nasser_ep_benchmarking_dataset:
  input: "results/ENCODE/ENCODE_Nasser2021_GRCh38.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_Nasser2021_GRCh38.tsv.gz"
  params:
     effect_size = "pctChange",
     min_pct_change = None,
     cell_type = "K562"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"

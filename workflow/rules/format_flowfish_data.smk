## Rules to download and format Engreitz lab Flow-FISH data for EP benchmarking

# download Engreitz lab Flow-FISH data
rule download_fulco:
  output: temp("results/ENCODE/FlowFISH_{sample}_data.tsv")
  params:
    url = lambda wildcards: config["download_urls"]["FlowFISH_" + wildcards.sample]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# reformat into EP benchmarking format
rule format_fulco:
  input: "results/ENCODE/FlowFISH_{sample}_data.tsv"
  output: "results/ENCODE/EPCrisprBenchmark_FlowFISH_{sample}_hg19.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/flowfish_dataset/format_flowfish.R"

# lift enhancer coordinates from hg19 to hg38 using UCSC's liftOver software    
rule liftover_flowfish_enhancers:
  input:
    results = "results/ENCODE/EPCrisprBenchmark_FlowFISH_{sample}_hg19.tsv.gz",
    chain = "results/ENCODE/liftover/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/ENCODE/liftover/FlowFISH_{sample}/enh_hg19.bed",
    hg38 = "results/ENCODE/liftover/FlowFISH_{sample}/enh_hg38.bed",
    unlifted = "results/ENCODE/liftover/FlowFISH_{sample}/enh_unlifted.bed"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "zcat {input.results} | "
    """awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $1":"$2"-"$3, 0, "."}}' | """
    "sort | uniq > {output.hg19} ; "
    "liftOver {output.hg19} {input.chain} {output.hg38} {output.unlifted}"

# liftover Gasperini2019 data from hg19 to hg38
rule liftover_flowfish_dataset:
  input:
    results = "results/ENCODE/EPCrisprBenchmark_FlowFISH_{sample}_hg19.tsv.gz",
    enh_hg38 = "results/ENCODE/liftover/FlowFISH_{sample}/enh_hg38.bed",
    annot_hg38 = "resources/gencode.v29.annotation.gtf.gz"
  output: "results/ENCODE/EPCrisprBenchmark_FlowFISH_{sample}_GRCh38.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/liftover_crispr_dataset.R"

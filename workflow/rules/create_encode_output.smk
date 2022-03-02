## Create input files for ENCODE distal regulation EP benchmarking pipelines

ruleorder: liftover_crispr_dataset > create_ep_benchmarking_dataset

# compile output files in ENCODE format
rule create_encode_dataset:
  input:
    results = "results/{sample}/output_{sd}gStd_{method}_{strategy}.tsv.gz",
    annot = "resources/gencode.v26lift37.annotation.gtf.gz",
    guide_targets = "resources/{sample}/guide_targets.tsv"
  output: "results/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}.tsv.gz"
  params:
    tss_min_dist = 1000,
    gene_ids = lambda wildcards: config["metadata"][wildcards.sample]["gene_ids"],
    tss_ctrl_tag = lambda wildcards: config["metadata"][wildcards.sample]["tss_ctrl_tag"],
    padj_threshold = config["diff_expr"]["padj_threshold"],
    cell_type = lambda wildcards: config["metadata"][wildcards.sample]["cell_type"],
    reference = lambda wildcards: config["metadata"][wildcards.sample]["reference"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_encode_dataset.R"
    
# extract minimal set of columns for EP benchmarking from full per element dataset    
rule create_ep_benchmarking_dataset:
  input: "results/ENCODE/ENCODE_{sample}_{sd}gStd_MAST_perCRE.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark_{sample}_{sd}gStd_unfiltered_{genome}.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  params:
     tss_to_dist = config["ep_benchmarking"]["dist_to_TSS"]
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"
    
# filter EP benchmarking datasets for minimum power
rule power_filter:
  input: "results/ENCODE/EPCrisprBenchmark_{sample}_{sd}gStd_unfiltered_{genome}.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark_{sample}_{sd}gStd_{pwr}pwrAt{es}effect_{genome}.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  params:
     remove_filtered_pairs = True
  script:
    "../scripts/encode_datasets/filter_ep_benchmarking_dataset.R"
    
# liftover CRISPRi datasets ------------------------------------------------------------------------

# download UCSC hg19 to hg38 liftover chain file
rule download_chain_file:
  output: "results/ENCODE/liftover/hg19ToHg38.over.chain.gz"
  params:
    url = config["download_urls"]["liftover_chain"]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# lift enhancer coordinates from hg19 to hg38 using UCSC's liftOver software    
rule liftover_enhancers:
  input:
    results = "results/ENCODE/EPCrisprBenchmark_Gasperini2019_{sd}gStd_unfiltered_hg19.tsv.gz",
    chain = "results/ENCODE/liftover/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/ENCODE/liftover/Gasperini2019_{sd}gStd_unfiltered/enh_hg19.bed",
    hg38 = "results/ENCODE/liftover/Gasperini2019_{sd}gStd_unfiltered/enh_hg38.bed",
    unlifted = "results/ENCODE/liftover/Gasperini2019_{sd}gStd_unfiltered/enh_unlifted.bed"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "zcat {input.results} | "
    """awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $1":"$2"-"$3, 0, "."}}' | """
    "sort | uniq > {output.hg19} ; "
    "liftOver {output.hg19} {input.chain} {output.hg38} {output.unlifted}"

# liftover Gasperini2019 data from hg19 to hg38
rule liftover_crispr_dataset:
  input:
    results = "results/ENCODE/EPCrisprBenchmark_Gasperini2019_{sd}gStd_unfiltered_hg19.tsv.gz",
    enh_hg38 = "results/ENCODE/liftover/Gasperini2019_{sd}gStd_unfiltered/enh_hg38.bed",
    annot_hg38 = "resources/gencode.v26.annotation.gtf.gz"
  output: "results/ENCODE/EPCrisprBenchmark_Gasperini2019_{sd}gStd_unfiltered_GRCh38.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark_Gasperini2019_{sd}gStd_{pwr}pwrAt{es}effect_GRCh38.tsv.gz",
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/liftover_crispr_dataset.R"

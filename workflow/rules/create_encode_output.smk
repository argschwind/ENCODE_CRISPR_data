## Create input files for ENCODE and distal regulation CRISPR benchmarking pipeline

ruleorder: liftover_crispr_dataset > create_encode_dataset

# function to get samples that require liftover from hg19 to GRCh38
def liftover_samples(config):
  genome_builds = config["encode_datasets"]["genome_build"].items()
  liftover_samples = list(dict(filter(lambda x: x[1] == "hg19", genome_builds)).keys())
  return(liftover_samples)

# compile output files in ENCODE format
rule create_encode_dataset:
  input:
    results = "results/{sample}/output_{sd}gStd_{method}_{strategy}.tsv.gz",
    annot = lambda wildcards: config["encode_datasets"]["annot"][wildcards.sample],
    guide_targets = "resources/{sample}/guide_targets.tsv"
  output: "results/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{genome}.tsv.gz"
  params:
    ignore_txs = lambda wildcards: config["encode_datasets"]["ignore_transcripts"][wildcards.sample],
    tss_min_dist = config["encode_datasets"]["dist_to_TSS"][0],
    gene_ids = lambda wildcards: config["metadata"][wildcards.sample]["gene_ids"],
    tss_ctrl_tag = lambda wildcards: config["metadata"][wildcards.sample]["tss_ctrl_tag"],
    padj_threshold = config["diff_expr"]["padj_threshold"],
    reference = lambda wildcards: config["metadata"][wildcards.sample]["reference"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_encode_dataset.R"
    
## Liftover CRISPRi datasets -----------------------------------------------------------------------

# lift enhancer coordinates from hg19 to hg38 using UCSC's liftOver software    
rule liftover_enhancers:
  input:
    results = "results/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz",
    chain = "resources/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg19.bed",
    hg38 = "results/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg38.bed",
    unlifted = "results/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_unlifted.bed"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "zcat {input.results} | "
    """awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $1":"$2"-"$3, 0, "."}}' | """
    "sort | uniq > {output.hg19} ; "
    "liftOver {output.hg19} {input.chain} {output.hg38} {output.unlifted}"
    
# liftover EP benchmarking dataset from hg19 to hg38
rule liftover_crispr_dataset:
  input:
    results = "results/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz",
    enh_hg38 = "results/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg38.bed",
    annot_hg38 = "resources/gencode.v26.annotation.gtf.gz"
  output: "results/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_GRCh38.tsv.gz"
  wildcard_constraints:
    sample = "|".join(liftover_samples(config))
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/liftover_crispr_dataset.R"    

## Create EPBenchmarking CRISPR data files ---------------------------------------------------------

# filter EP benchmarking datasets for distance to TSS and minimum power
rule filter_crispr_dataset:
  input: "results/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{genome}.tsv.gz"
  output: temp("results/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{pwr}pwrAt{es}effect_{genome}.tsv.gz")
  params:
    tss_to_dist = config["encode_datasets"]["dist_to_TSS"],
    remove_filtered_pairs = False
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/filter_crispr_dataset.R"

# convert ENCODE format files to EPBenchmarking format files
rule create_ep_benchmarking_dataset:
  input: "results/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_{sd}gStd_MAST_perCRE_{pwr}pwrAt{es}effect_{genome}.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_{sample}_{sd}gStd_{pwr}pwrAt{es}effect_{genome}.tsv.gz"
  params:
    effect_size = "logFC",
    min_pct_change = None,
    cell_type = lambda wildcards: config["metadata"][wildcards.sample]["cell_type"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"    

## Create ensemble dataset -------------------------------------------------------------------------

# create ensembl CRISPR dataset in both ENCODE format    
rule create_ensemble_encode:
  input:
    Fulco2019 = "results/ENCODE/ENCODE_Fulco2019_GRCh38.tsv.gz",
    Gasperini2019 = "results/ENCODE/EPCrisprBenchmark/ENCODE_Gasperini2019_0.13gStd_MAST_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz",
    Schraivogel2020 = "results/ENCODE/EPCrisprBenchmark/ENCODE_TAPseq_0.13gStd_MAST_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz"
  output: "results/ENCODE/ENCODE_CombinedData_GRCh38.tsv.gz"
  params:
    effect_size = {"Fulco2019":"pctChange", "Gasperini2019":"logFC", "Schraivogel2020":"logFC"}
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ensemble_dataset.R"

# convert ensembl CRISPR dataset from ENCODE to EPBenchmarking format file  
rule create_ensemble_epbenchmarking:
  input: "results/ENCODE/ENCODE_CombinedData_GRCh38.tsv.gz"
  output: "results/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_CombinedData_GRCh38.tsv.gz"
  params:
    effect_size = "pctChange",
    min_pct_change = None,
    cell_type = "K562"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/encode_datasets/create_ep_benchmarking_dataset.R"   

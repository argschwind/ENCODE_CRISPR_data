
rule download_processed_hct116:
  output: "results/HCT116/hct116_hg19.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/hct116_dataset/download_hct116_data.R"
  
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

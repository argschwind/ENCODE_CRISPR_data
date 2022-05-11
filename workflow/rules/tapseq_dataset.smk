## Rules to create differential expression input data for Schraivogel et al TAP-seq data

# Create input data for TAP-seq experiments --------------------------------------------------------

# match guide sequences to genome sequence to get binding sites and create hg19 guide targets file
rule create_tapseq_guide_targets:
  input:
    vectors = "resources/TAPseqChr{chr}/cropseq_vectors_chr{chr}_screen.fasta",
    annot = "resources/gencode.v29.annotation.gtf.gz"
  output:
    guide_targets = "resources/TAPseqChr{chr}/guide_targets_hg19.tsv",
    guides_bed = "resources/TAPseqChr{chr}/guide_positions_hg19.bed",
    targets_bed = "resources/TAPseqChr{chr}/guide_targets_hg19.bed"
  params:
    bsgenome = "BSgenome.Hsapiens.UCSC.hg19",
    upstream_seq = "GGCTTTATATATCTTGTGGAAAGGACGAAACACCG",
    downstream_seq = "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTT"
  threads: 10
  conda: "../envs/r_create_tapseq_input.yml"
  script:
    "../scripts/tapseq_dataset/map_tapseq_guides.R"
  
# liftover guide coordinates from hg19 to hg38
rule liftover_guides:
  input:
    guides_hg19 = "resources/{sample}/guide_positions_hg19.bed",
    chain = "resources/hg19ToHg38.over.chain.gz"
  output:
    guides_hg38 = "resources/{sample}/guide_positions_hg38.bed",
    unlifted = "resources/{sample}/guide_positions_hg19_unlifted.bed"
  conda: "../envs/r_create_tapseq_input.yml"
  shell:
    "liftOver {input.guides_hg19} {input.chain} {output.guides_hg38} {output.unlifted}"

# liftover target coordinates from hg19 to hg38
rule liftover_targets:
  input:
    targets_hg19 = "resources/{sample}/guide_targets_hg19.bed",
    chain = "resources/hg19ToHg38.over.chain.gz"
  output:
    targets_hg38 = "resources/{sample}/guide_targets_hg38.bed",
    unlifted = "resources/{sample}/guide_targets_hg19_unlifted.bed"
  conda: "../envs/r_create_tapseq_input.yml"
  shell:
    "liftOver {input.targets_hg19} {input.chain} {output.targets_hg38} {output.unlifted}"

# liftover guide targets file
rule liftover_guide_targets_file:
  input:
    guide_targets_hg19 = "resources/{sample}/guide_targets_hg19.tsv",
    guides_hg38 = "resources/{sample}/guide_positions_hg38.bed",
    targets_hg38 = "resources/{sample}/guide_targets_hg38.bed"
  output: "resources/{sample}/guide_targets_hg38.tsv"
  conda: "../envs/r_create_tapseq_input.yml"
  script:
    "../scripts/tapseq_dataset/liftover_guide_targets.R"

# create Perturb-seq SCE objects for TAP-seq samples
rule create_sce_tapseq:
  input:
    dge = "resources/{sample}/dge.txt",
    pert_status = "resources/{sample}/perturb_status.txt",
    annot = "resources/gencode.v29.annotation.gtf.gz",
    guide_targets = "resources/{sample}/guide_targets_hg38.tsv"
  output: temp("resources/{sample}/perturb_sce_unfilt.rds")
  params:
    vector_pattern = "^CROPseq_dCas9_DS_.+$"
  log: "resources/{sample}/logs/create_sce_tapseq.log"
  conda: "../envs/r_create_tapseq_input.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/tapseq_dataset/create_sce_tapseq.R"
    
# add 10x lanes to SCE colData and filter out cells from bad lanes
rule filter_tapseq_10x_lanes:
  input: "resources/{sample}/perturb_sce_unfilt.rds"
  output: "resources/{sample}/perturb_sce.rds"
  wildcard_constraints:
    sample = "TAPseqChr.+"
  params:
    remove_lanes = lambda wildcards: config["process_tapseq"]["remove_lanes"][wildcards.sample]
  log: "resources/{sample}/logs/filter_tapseq_10x_lanes.log"
  conda: "../envs/r_create_tapseq_input.yml"
  script:
    "../scripts/tapseq_dataset/filter_tapseq_10x_lanes.R"
    
# Combine output from TAP-seq experiments ----------------------------------------------------------

# combine guide targets file from both experiments
rule combine_tapseq_guide_targets:
  input:
    chr8 = "resources/TAPseqChr8/guide_targets_hg38.tsv",
    chr11 = "resources/TAPseqChr11/guide_targets_hg38.tsv"
  output:
    guide_targets = "resources/TAPseq/guide_targets.tsv",
    guides_bed = "resources/TAPseq/guide_positions.bed",
    targets_bed = "resources/TAPseq/guide_targets.bed"
  conda: "../envs/r_create_tapseq_input.yml"
  script:
    "../scripts/tapseq_dataset/combine_tapseq_guide_targets.R"

# combine TAP-seq differential expression results from both experiments
rule combine_tapseq_de_results:
  input: 
    chr8 = "results/TAPseqChr8/diff_expr/output_{method}_{strategy}.tsv.gz",
    chr11 = "results/TAPseqChr11/diff_expr/output_{method}_{strategy}.tsv.gz"
  output: "results/TAPseq/diff_expr/output_{method}_{strategy}.tsv.gz"
  params: 
    p_adj_method = config["diff_expr"]["padj_method"]
  conda: "../envs/r_create_tapseq_input.yml"
  script:
    "../scripts/tapseq_dataset/combine_de_results.R"
    
# compute power across both TAP-seq experiments from simulation results
rule compute_tapseq_power:
  input:
    chr8 = expand("results/TAPseqChr8/power_sim/rep{rep}_output_{{effect}}_{{sd}}gStd_{{method}}_{{strategy}}.tsv.gz",
      rep = range(1, config["power_simulations"]["rep"] + 1)),
    chr11 = expand("results/TAPseqChr11/power_sim/rep{rep}_output_{{effect}}_{{sd}}gStd_{{method}}_{{strategy}}.tsv.gz",
      rep = range(1, config["power_simulations"]["rep"] + 1))
  output:
    "results/TAPseq/power_sim/power_{effect}_{sd}gStd_{method}_{strategy}.tsv.gz"
  params:
    p_adj_method = config["diff_expr"]["padj_method"],
    pval_threshold = config["diff_expr"]["padj_threshold"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/tapseq_dataset/compute_tapseq_power.R"

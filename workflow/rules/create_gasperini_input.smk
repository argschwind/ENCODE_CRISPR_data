## Rules to create differential expression input data for Gasperini et al Perturb-seq data

# download Gasperini et al data --------------------------------------------------------------------

# download main Gasperini et al. 2019 data object
rule download_gasperini_cds:
  output: temp("resources/Gasperini2019/GSE120861_at_scale_screen.cds.rds")
  params:
    url = config["download_urls"]["gasperini_cds"]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output}.gz {params.url} ; gunzip {output}.gz"
    
# download annotations used by Gasperini et al.
rule download_gencode_annotations:
  output: "resources/gencode.v26lift37.annotation.gtf.gz"
  params:
    url = config["download_urls"]["gencode_v26lift37"]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# download and process DNase-seq data --------------------------------------------------------------
    
# download hg19 DNase-seq file also used for ABC predictions
rule download_dnase:
  output:
    bam = "resources/DNase/{dnase_bam}/{dnase_bam}.bam",
    bai = "resources/DNase/{dnase_bam}/{dnase_bam}.bam.bai"
  params:
    url = lambda wildcards: config["download_urls"][wildcards.dnase_bam]
  conda: "../envs/r_create_gasperini_input.yml"
  shell:
    "wget -O {output.bam} {params.url} ; samtools index {output.bam}"

# perform peak calling using MACS2
rule call_peaks:
  input:
    bam = "resources/DNase/{dnase_bam}/{dnase_bam}.bam",
    bai = "resources/DNase/{dnase_bam}/{dnase_bam}.bam.bai"
  output: 
    peaks = "resources/DNase/{dnase_bam}/{dnase_bam}_peaks.narrowPeak",
    summits = "resources/DNase/{dnase_bam}/{dnase_bam}_summits.bed"
  log: "resources/DNase/{dnase_bam}/logs/call_peaks.log"
  params:
    outdir = "resources/DNase/{dnase_bam}",
  conda: "../envs/r_create_gasperini_input.yml"
  shell:
    "macs2 callpeak "
    "-t {input.bam} "
    "-n {wildcards.dnase_bam} "
    "-f BAM "
    "-g hs "
    "-p .1 "
    "--call-summits "
    "--outdir {params.outdir} "
    "2> {log}"
    
# create DNase-seq summit regions by extending peaks by a specified amount
rule create_summit_regions:
  input: 
    summits = "resources/DNase/{dnase_bam}/{dnase_bam}_summits.bed",
    chr_sizes = "resources/hg19_chr_sizes"
  output: "resources/DNase/{dnase_bam}/{dnase_bam}_summit_regions.bed"
  params:
    extend_bp = 150
  conda: "../envs/r_create_gasperini_input.yml"
  shell:
    "bedtools slop -b {params.extend_bp} -i {input.summits} -g {input.chr_sizes} | "
    """awk '{{print $0"\\t*"}}' > {output}"""
   
# create Gasperini differential expression analysis input data -------------------------------------

# reformat remapped guides to bed file with only remapped coordinates (exclude non-targeting ctrls)
rule create_guides_bed:
  input: 
    guides = "resources/Gasperini2019/GasperiniGuidePositions.tsv",
    chr_sizes = "resources/hg19_chr_sizes"
  output: "resources/Gasperini2019/guide_positions.bed"
  conda: "../envs/r_create_gasperini_input.yml"
  shell:
    """awk 'BEGIN {{FS = "\\t"}} ; NR != 1 && $10 != "NA" && $7 != "NTC" """
    """{{print $9"\\t"$10"\\t"$11"\\t"$1"\\t"$12"\\t"$13"\\t"$2}}' {input.guides} | """
    "bedtools sort -i stdin -g {input.chr_sizes} > {output}"

# create gRNA targets file from gRNA and perturbation target coordinates
rule create_targets_file:
  input:
    grnas = "resources/Gasperini2019/guide_positions.bed",
    targets = "resources/DNase/ENCFF030DCL/ENCFF030DCL_summit_regions.bed"
  output:
    grna_targets_file = "resources/Gasperini2019/guide_targets.tsv",
    n_overlaps = "resources/Gasperini2019/guide_targets_overlaps.tsv"
  log: "resources/Gasperini2019/logs/create_targets_file.log"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_grna_targets_file.R"
    
# convert to summarized experiment object like for TAP-seq data, with added gene annotation
rule create_sce_gasperini:
  input:
    cds = "resources/Gasperini2019/GSE120861_at_scale_screen.cds.rds",
    annot = "resources/gencode.v26lift37.annotation.gtf.gz",
    guide_targets = "resources/Gasperini2019/guide_targets.tsv"
  output: "resources/Gasperini2019/perturb_sce.rds"
  log: "resources/Gasperini2019/logs/create_sce_gasperini.log"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/gasperini_dataset/create_sce_gasperini.R"

# alternative guide targets ------------------------------------------------------------------------

# download Gasperini et al. 2019 gRNA library
rule download_gasperini_grnas:
  output: "resources/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx"
  params:
    url = "https://ars.els-cdn.com/content/image/1-s2.0-S009286741831554X-mmc2.xlsx"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# create gRNA targets file based on target definitions from Gasperini et al. 2019
rule create_gasperini_guide_targets:
  input:
    guides = "resources/Gasperini2019/guide_positions.bed",
    targets = "resources/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx"
  output:
    guide_targets = "resources/Gasperini2019/gasperini_targets/guide_targets.tsv",
    targets_bed = "resources/Gasperini2019/gasperini_targets/gasperini_targets.bed"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/gasperini_dataset/create_gasperini_guide_targets.R"
    
# create guide targets file based on Joe's approach of merging proximal guides into elements
rule create_joe_guide_targets:
  input: "resources/Gasperini2019/guide_positions.bed"
  output: 
    guide_targets = "resources/Gasperini2019/joe_targets/guide_targets.tsv",
    targets_bed = "resources/Gasperini2019/joe_targets/joe_targets.bed"
  params:
    extend_bp = 75
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/gasperini_dataset/create_joe_guide_targets.R"

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
    
# download Gasperini et al. 2019 gRNA library
rule download_gasperini_grnas:
  output: "resources/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx"
  params:
    url = "https://ars.els-cdn.com/content/image/1-s2.0-S009286741831554X-mmc2.xlsx"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"
    
# download gencode annotations
rule download_gencode_annotations:
  output: "resources/{annot}.annotation.gtf.gz"
  params:
    url = lambda wildcards: config["download_urls"][wildcards.annot]
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

# perform peak calling using MACS2 (high p-value cutoff leads to many weak peaks)
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
    
# quantify DNase-seq seq reads in peaks
rule quantify_dnase:
  input:
    peaks = "resources/DNase/{dnase_bam}/{dnase_bam}_peaks.narrowPeak",
    bam = "resources/DNase/{dnase_bam}/{dnase_bam}.bam",
    chrs = "resources/hg19_chr_sizes"
  output: "resources/DNase/{dnase_bam}/{dnase_bam}_peaks.counts"
  log: "resources/DNase/{dnase_bam}/logs/quantify_dnase.log"
  conda: "../envs/r_create_gasperini_input.yml" 
  shell:
    "bedtools sort -faidx {input.chrs} -i {input.peaks} | "
    "bedtools coverage -sorted -a stdin -b {input.bam} -g {input.chrs} -counts > {output}"
    
# create candidate CREs by selecting top peaks, extract and extend summits, and merge overlaps
rule create_candidate_cres:
  input: 
    peaks = "resources/DNase/{dnase_bam}/{dnase_bam}_peaks.counts",
    chrs = "resources/hg19_chr_sizes"
  output: "resources/DNase/{dnase_bam}/{dnase_bam}_candidate_cres.bed"
  params:
    top_n = 150000,
    peak_extend = 250
  conda: "../envs/r_create_gasperini_input.yml"
  shell:
    "bedtools sort -i {input.peaks} -faidx {input.chrs} | "
    "bedtools merge -i stdin -c 11 -o max | "
    "sort -nr -k 4 | awk '(NR <= 150000)' | "
    "bedtools intersect -b stdin -a {input.peaks} -wa | "
    """awk 'BEGIN {{OFS="\\t"}} {{print $1, $2 + $10, $2 + $10}}' | """ 
    "bedtools slop -i stdin -b {params.peak_extend} -g {input.chrs} | "
    "bedtools sort -i stdin -faidx {input.chrs} | "
    "bedtools merge -i stdin | "
    """awk 'BEGIN {{OFS="\\t"}} {{print $1, $2, $3, $1":"$2"-"$3, "0", "."}}' > {output}"""

   
# create Gasperini differential expression analysis input data -------------------------------------

# reformat remapped guides to bed file with only remapped coordinates (exclude non-targeting ctrls)
rule create_guides_bed:
  input: 
    guides = "resources/Gasperini2019/GasperiniGuidePositions.tsv",
    chrs = "resources/hg19_chr_sizes"
  output: "resources/Gasperini2019/guide_positions.bed"
  conda: "../envs/r_create_gasperini_input.yml"
  script:
    "../scripts/gasperini_dataset/create_gasperini_guides_bed.R"

# extract Gasperini positive control targets and guides
rule gasperini_pos_ctrls:
  input:
    guides = "resources/Gasperini2019/guide_positions.bed",
    targets = "resources/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx"
  output: temp("resources/Gasperini2019/gasperini_pos_ctrls.tsv")
  conda: "../envs/r_create_gasperini_input.yml"
  script:
    "../scripts/gasperini_dataset/gasperini_pos_ctrls.R"
    
# create guide targets file by overlapping gRNAs with candidate CRE coordinates from DNase peaks
rule create_guide_targets_file:
  input:
    guides = "resources/Gasperini2019/guide_positions.bed",
    targets = "resources/DNase/ENCFF030DCL/ENCFF030DCL_candidate_cres.bed",
    known_targets = "resources/Gasperini2019/gasperini_pos_ctrls.tsv"
  output:
    guide_targets_file = "resources/Gasperini2019/guide_targets.tsv",
    targets_bed = "resources/Gasperini2019/guide_targets.bed",
    targets_per_guide_file = "resources/Gasperini2019/guide_targets_overlaps.tsv"
  log: "resources/Gasperini2019/logs/create_targets_file.log"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_guide_targets_file.R"
    
# convert to summarized experiment object like for TAP-seq data, with added gene annotation
rule create_sce_gasperini:
  input:
    cds = "resources/Gasperini2019/GSE120861_at_scale_screen.cds.rds",
    annot = "resources/gencode.v26lift37.annotation.gtf.gz",
    guide_targets = "resources/Gasperini2019/guide_targets.tsv"
  output: "resources/Gasperini2019/perturb_sce.rds"
  log: "resources/Gasperini2019/logs/create_sce_gasperini.log"
  conda: "../envs/r_create_gasperini_input.yml"
  resources:
    mem = "32G"
  script:
    "../scripts/gasperini_dataset/create_sce_gasperini.R"

# alternative guide targets ------------------------------------------------------------------------

# create gRNA targets file based on target definitions from Gasperini et al. 2019
rule create_gasperini_guide_targets:
  input:
    guides = "resources/Gasperini2019/guide_positions.bed",
    targets = "resources/Gasperini2019/1-s2.0-S009286741831554X-mmc2.xlsx"
  output:
    guide_targets = "resources/Gasperini2019/gasperini_targets/guide_targets.tsv",
    targets_bed = "resources/Gasperini2019/gasperini_targets/gasperini_targets.bed"
  conda: "../envs/r_create_gasperini_input.yml"
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
  conda: "../envs/r_create_gasperini_input.yml"
  script:
    "../scripts/gasperini_dataset/create_joe_guide_targets.R"

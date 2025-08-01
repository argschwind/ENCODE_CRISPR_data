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
  params:
    targets_type = "enh"
  log: "resources/Gasperini2019/logs/create_targets_file.log"
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_guide_targets_file.R"
    
# create Perturb-seq SCE object for Gasperini dataset
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

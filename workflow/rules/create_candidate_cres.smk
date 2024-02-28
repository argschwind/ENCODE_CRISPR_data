## Rules to download DNase-seq data and create candidate CREs
    
# download DNase-seq bam files
rule download_dnase:
  output:
    bam = temp("resources/DNase/{dnase_bam}/{dnase_bam}.bam"),
    bai = temp("resources/DNase/{dnase_bam}/{dnase_bam}.bam.bai")
  params:
    url = lambda wildcards: config["dnase_bam"][wildcards.dnase_bam]["url"]
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
    chrs = lambda wildcards: config["dnase_bam"][wildcards.dnase_bam]["chrs"]
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
    chrs = lambda wildcards: config["dnase_bam"][wildcards.dnase_bam]["chrs"]
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

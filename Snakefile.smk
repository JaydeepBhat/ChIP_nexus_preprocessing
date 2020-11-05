"""
Author: Jaydeep Bhat
Affiliation: TUM school of Life Sciences Weihenstephan, Freising; Date: March 2020; Run: snakemake
------------------
Modified from original author: Melanie Weilert
Affiliation: Stowers Institute; Aim: ChIP-nexus pipeline (lab code); Date: September 2019; Run: snakemake
------------------
Modifications (fixes) to the original code:
- rule bam_to_granges (syntax error)
- rule bam_to_granges (shell command error)
- rule bam_to_granges (shell command addition)
- rule generate_bed_file
- rule generate_clean_bed_file
------------------
Environment: (conda): GalwayGirlInTheSummerOf69
------------------
Original requirements:
------------------
- R == 3.6.1
- Bioconductor >= 3.2
- bowtie == 1.2.3
- gzcat
- cutadapt == 2.6
- samtools == 1.9

Main target rules:
------------------
- fastq_preprocess_single_end: trim and record barcoding (fixed and unique) from ChIP-nexus single end sequences
- fastq_trim_adapters: trim ChIP-nexus adapter sequences from readsIndices
- fastq_align: align FASTQ reads to the specified genome
# - bam_deduplicate: deduplicate BAM files
- bam_to_granges: convert deduplicated BAM files to GRanges RDS files, save deduplicated BAM as intermediate file
- granges_to_bw: convert ChIP-nexus GRanges files to bigwig files
- bw_normalize: normalize ChIP-nexus bw files to RPM
- rule generate_clean_bed_file: convert peaks into bed format for downstream analysis (e.g. motifs)
"""

#Setup
import csv
import os
from itertools import product

##############################################################################
#Required inputs to run Snakemake...raw FASTQ files should follow {factor}_nexus_{rep} naming convention.
##############################################################################

FACTORS=["sample"] #factors to process
REPS=[1] #replicates for each factor
FIXED_BARCODES='CTGA,TGAC,GACT,ACTG,CAGT,TCAG,GTCA,AGTC' #rev comp of oligo fixed barcodes in ChIP-nexus experiment
BOWTIE_INDEXES='bowtie_index/BASENAME' # change "BASENAME" in this string to pre-build Bowtie index basename (read more info in the file "readme.info.txt")

##############################################################################
#If you have missing replicates, edit "forbidden" variable to add missing reps.
##############################################################################
def filter_combinator(combinator, blacklist):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            # Use frozenset instead of tuple
            # in order to accomodate
            # unpredictable wildcard order
            if frozenset(wc_comb) not in blacklist:
                yield wc_comb
    return filtered_combinator

forbidden = {
    # frozenset({("factor", "klf4"), ("rep", 3)}), #e.g
    }
filtered_samples = filter_combinator(product, forbidden)

##############################################################################
# Create rule to determine outputs based on FACTOR and REP above.
##############################################################################
rule all:
    input:
       expand("bw/{factor}_nexus_{rep}_filtered_positive.bw", factor=FACTORS, rep=REPS),
       expand("bw/{factor}_nexus_{rep}_filtered_normalized_positive.bw", factor=FACTORS, rep=REPS),
       expand("macs2/{factor}_nexus_{rep}_peaks.narrowPeak", factor=FACTORS, rep=REPS),
       expand("bed/{factor}_nexus_{rep}_clean.bed", factor=FACTORS, rep=REPS)

       # the commented out files should be used instead if you were missing replicates
       # expand("bw/{factor}_nexus_{rep}_filtered_positive.bw", filtered_samples, factor=FACTORS, rep=REPS),
       # expand("bw/{factor}_nexus_{rep}_filtered_normalized_positive.bw", filtered_samples, factor=FACTORS, rep=REPS),
       # expand("macs2/{factor}_nexus_{rep}_peaks.narrowPeak", filtered_samples, factor=FACTORS, rep=REPS),

##############################################################################
# Rules in pipeline are reported below
##############################################################################

# Trim and record ChIP-nexus barcodes (fixed and unique) from the FASTQ
rule fastq_preprocess_single_end:
    input:
        "/PATH_TO_YOUR_RAW_DATA_FASTQ_FILES/fastq/raw/{factor}_nexus_{rep}.fastq.gz",
    output:
        "fastq/processed/{factor}_nexus_{rep}_processed.fastq.gz",
    params:
        barcodes = FIXED_BARCODES,
        pretrim_length = 50, #trim sequence to this length
        minimum_fragment_length = 22, #keep fragments longer than this length
        unique_barcode_length = 5, #length of unique barcodes
        nworkers = 4, #numer of parallel processes
        chunksize = 1000, #how many reads to process at a time
    message:
        "Preprocessing single end ChIP-nexus sample..."
    shell:
        "Rscript scripts/preprocess_fastq.r -f {input} -t {params.pretrim_length} \
        -k {params.minimum_fragment_length} -b {params.barcodes} -r {params.unique_barcode_length} \
        -p {params.nworkers} -c {params.chunksize} -o {output}"

# Remove ChIP-nexus adapters from all reads
rule fastq_trim_adapters: #Cutadapt should be version 1.8.1
    input:
        "fastq/processed/{factor}_nexus_{rep}_processed.fastq.gz",
    output:
        "fastq/trimmed/{factor}_nexus_{rep}_processed_trimmed.fastq.gz",
    params:
        minimum_overlap = 4, #minimum overlap allowed
        minimum_fragment_length = 22, #keep fragments longer than this length
        maximum_error_rate = .2, #maximum error rate allowed for adapter searches
    message:
        "Preprocessing single end ChIP-nexus sample..."
    shell:
        "cutadapt -m {params.minimum_fragment_length} -O {params.minimum_overlap} \
        -e {params.maximum_error_rate} --quiet \
        -a AGATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC -o {output} {input}"

# Align ChIP-nexus samples using bowtie1
rule fastq_align:
    input:
        "fastq/trimmed/{factor}_nexus_{rep}_processed_trimmed.fastq.gz",
    output:
        "bam/{factor}_nexus_{rep}.bam",
    params:
        indexes = BOWTIE_INDEXES, #indexes are faster versions of .fa files, specific to bowtie1
        nworkers = 4, #numer of parallel processes
        alignments_reported = 1, #how many valid alignments reported
        multialignments_allowed = 1, #suppress all alignments if more than these reportable alignments exist
        mismatches_allowed = 2, #mismatches allowed in alignments
    message:
        "Aligning ChIP-nexus sample...",
    shell:
        "gzip -cd {input} | bowtie -S -p {params.nworkers} --chunkmbs 512 \
        -k {params.alignments_reported} -m {params.multialignments_allowed} \
        -v {params.mismatches_allowed} --best --strata \
        {params.indexes} - | samtools view -F 4 -Sbo {output} -"

#Convert BAM to rdata
rule bam_to_granges:
    input:
        "bam/{factor}_nexus_{rep}.bam",
    output:
        "rdata/{factor}_nexus_{rep}_filtered.granges.rds",
        "bam/{factor}_nexus_{rep}_filtered.bam",
    params:
        name = "rdata/{factor}_nexus_{rep}",
    message: "Converting BAM to GRanges..."
    run:
        shell("Rscript scripts/process_bam.r -f {input} -n {params.name} -o TRUE")
        shell("mv {params.name}*.bam bam -v")
        shell("mv {params.name}*.bai bam -v")
        shell("mv {params.name}.granges.rds {params.name}_filtered.granges.rds -v")

#Convert rdata to bw
rule granges_to_bw:
    input:
        "rdata/{factor}_nexus_{rep}_filtered.granges.rds",
    output:
        "bw/{factor}_nexus_{rep}_filtered_positive.bw",
    params:
        fac = "{factor}",
        r = "{rep}",
    message: "Converting GRanges to bw..."
    run:
        shell("Rscript scripts/split_granges_by_strand.r -r {input}")
        shell("mv rdata/{params.fac}_nexus_{params.r}_filtered_*.bw bw -v")

rule bw_normalize:
    input:
        "rdata/{factor}_nexus_{rep}_filtered.granges.rds",
    output:
        "bw/{factor}_nexus_{rep}_filtered_normalized_positive.bw",
    params:
        name = "bw/{factor}_nexus_{rep}_filtered_normalized",
    message: "Normalizing ChIP-nexus data..."
    shell:
        "Rscript scripts/generating_normalized_tracks_from_gr.r -f {input} -n {params.name}"

#Run MACS2 on all replicates
rule run_nexus_macs2_individual:
    input:
        "bam/{factor}_nexus_{rep}_filtered.bam"
    output:
        "macs2/{factor}_nexus_{rep}_peaks.narrowPeak"
    params:
        name = "{factor}_nexus_{rep}",
        outdir = "macs2"
    message: "Running MACS2 on individual BAM files..."
    shell:
        "mkdir tmp; macs2 callpeak -t {input} --tempdir tmp --qvalue 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all -g mm -f BAM --outdir {params.outdir} -n {params.name}"

# generate bed files from narrow peaks for motif analysis
rule generate_bed_file:
    input:
        "macs2/{factor}_nexus_{rep}_peaks.narrowPeak"
    output:
        temp("bed/{factor}_nexus_{rep}.bed")
    params:
        name = "{factor}_nexus_{rep}",
        outdir = "bed"
    message: "What is in the bed... Motifs"
    shell:
        "sort -k8,8nr {input} | cut -f 1,2,3 > {output}"

rule generate_clean_bed_file:
    input:
        "bed/{factor}_nexus_{rep}.bed"
    output:
        "bed/{factor}_nexus_{rep}_clean.bed"
    params:
        name = "{factor}_nexus_{rep}",
        outdir = "bed"
    message: "clean the bed... for Motifs"
    shell:
        "less -S {input} | grep chr > {output}"

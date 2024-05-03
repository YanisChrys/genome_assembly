"""
Snakemake Pipeline for HiFi Genome Assemblies
This pipeline processes HiFi Genome Assemblies
             
Version: 1.0
Author: Ioannis Chrysostomakis
Email: john.chrysost@gmail.com
Date Created: YYYY-MM-DD
Last Modified: YYYY-MM-DD

License: MIT

Dependencies: 
tested on snakemake 7.24.0
Acknowledgements: 

Notes: Make sure the configuration file is properly set up 
and make sure to read the README file before attempting to run this pipeline.
Don't panic and read the comments :)
"""

### Load Configurations ###

import re
import os
configfile: "config/config.yaml"

# Save input from config as variables for easy access
#ENZYMES=config["enzymes"]
SUBREADS = config["subreads"]
HIFI = config["hifi"]
HIC1 = config["hic1"]
HIC2 = config["hic2"]
STOP_AT_HIFIASM = config["STOP_AT_HIFIASM"]
DO_STATS = config["DO_STATS"]
DO_BUSCO = config["DO_BUSCO"]
DO_HIC = config["DO_HIC"]
DO_CCS = config["DO_CCS"]
DO_DEDUP = config["DO_DEDUP"]
DO_ADAPT_FILT = config["DO_ADAPT_FILT"]
INCLUDE_HIC = config["INCLUDE_HIC"]
CHUNK_NMB = list(range(1,config["ccs"]["chunks"]+1))

### Constants ###

# remove special characters and dots from prefix but keep underscores
PREFIX = re.sub('[^a-zA-Z0-9 \_]', '', config["file_prefix"])

############ Libraries of assemblies ############

# snakemake is good at running on simiarly named files that sit on the same folder.
# Here, we have files on different folders with different names and snakemake is 
# less helpful when handling those. We bypass that by forcing it to do a 
# for-loop over all the files and using the basename of each file as an identifier.

# Additionally, we need to run different programs on different programs.
# For example, purge_dups needs to run on the hifiasm output but assembly statistics
# need to run on all fasta files. So we need different libraries for each task.

# The same is true for all the different iterations of our reads
# (pre- and post- decontamination and correction)

HIFIASM_OPTIONS = config["HIFIASM_OPTIONS"]

HIFIASM_PURGE = config["HIFIASM_PURGE"]

# create filenames with all hifiasm option combinations
HIFIASM_FILES = [
    "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_" + ops + "_" + l + "_primary.fasta" 
    for ops in HIFIASM_OPTIONS.keys() for l in HIFIASM_PURGE
] + [
    "RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_" + ops + "_" + l + "_alternate.fasta" 
    for ops in HIFIASM_OPTIONS.keys() for l in HIFIASM_PURGE
]
HIFIASM_BASENAMES=[os.path.basename(f).split('.')[0] for f in HIFIASM_FILES]

HIFIASM_MERQURY_BASENAMES = list(set(
    [os.path.basename(f).rsplit('_', 1)[0] for f in HIFIASM_FILES]
))

PURGE_DUP_FILES = expand("RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta", basename=HIFIASM_BASENAMES)
PURGE_DUP_BASENAMES = [os.path.basename(f).split('.')[0] for f in PURGE_DUP_FILES]

# yahs (or future: other) scaffoldsed assemblies
SCAF_FILES = [
    "RESULTS/HIC/YAHS/" + PREFIX + "_" + ops + "_" + l + "_primary_yahs_scaffolds_final.fa" 
    for ops in HIFIASM_OPTIONS.keys() for l in HIFIASM_PURGE
] + [
    "RESULTS/HIC/YAHS/" + PREFIX + "_" + ops + "_" + l + "_alternate_yahs_scaffolds_final.fa" 
    for ops in HIFIASM_OPTIONS.keys() for l in HIFIASM_PURGE
]
SCAF_BASENAMES=[os.path.basename(f).split('.')[0] for f in SCAF_FILES]

# all assemblies
ASSEMBLY_FILES=HIFIASM_FILES

if DO_DEDUP:
    ASSEMBLY_FILES=ASSEMBLY_FILES + PURGE_DUP_FILES
if DO_HIC:
    ASSEMBLY_FILES=ASSEMBLY_FILES + SCAF_FILES
ASSEMBLY_BASENAMES = [os.path.basename(f).split('.')[0] for f in ASSEMBLY_FILES]



## Do the same for the read files

# reads will be decontaminated regardless of whether they are produced 
HIFI_READ_FILES=["DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz"]
if HIFI: 
    HIFI_READ_FILES = HIFI_READ_FILES + [HIFI]

# if we give subreads as input, we still want to run stats on the hifi that get created during the run
if SUBREADS: 
    HIFI_READ_FILES = HIFI_READ_FILES + ["DATA/" + PREFIX + ".fastq.gz"]

HIC_READ_FILES=[]
if HIC1 and HIC2: 
    HIC_READ_FILES = HIC_READ_FILES + [HIC1] + [HIC2]
if not HIC1 and not HIC2: 
    HIC_READ_FILES = ["empty","empty"]

# Extracting basenames from HIFI_READ_FILES & HIC_READ_FILES
# these will be used to identify the different statistics used on each file
HIFI_BASENAMES = [os.path.basename(f).split('.')[0] for f in HIFI_READ_FILES]
HIC_BASENAMES = [os.path.basename(f).split('.')[0] for f in HIC_READ_FILES]

hic_pattern = "|".join(HIC_BASENAMES)

wildcard_constraints:
    chunknumber="[0-9]+",
    hic=hic_pattern,
    basename=".*_.*"



###  rule all input ###

# Snakemake requires you provide it with all the "terminal"
# output of your pipeline.
# This pipeline allows for the deactivation or activation
# of different parts of it. This means that not all rules can run 
# every time, ie not every run will have the same terminal files. 
# Below we conditionally add files to the rule all 
# depending on whether their rules are expected to run.
# Takes advantage of the fact the snakemake input is a list of strings.

OUTPUT_FILES = []
 
if STOP_AT_HIFIASM:
    DO_STATS = False
    DO_BUSCO = False
    DO_HIC = False
    DO_DEDUP = False
    DO_ADAPT_FILT = False
    hifiasm = expand("RESULTS/GENOME_ASSEMBLY/{basename}.fasta", basename=HIFIASM_BASENAMES)
    OUTPUT_FILES = OUTPUT_FILES + hifiasm
if DO_STATS:
    #assembly stats for all assemblies
    # read stats
    katgc1 = "RESULTS/STATISTICS/FASTK/" + PREFIX + ".st.png"
    katgc2 = "RESULTS/STATISTICS/FASTK/" + PREFIX + ".fi.png"
    katgc3 = "RESULTS/STATISTICS/FASTK/" + PREFIX + ".ln.png"
    smudgeplot = "RESULTS/STATISTICS/FASTK/" + PREFIX + ".png"
    genescope1 = "RESULTS/STATISTICS/FASTK/" + PREFIX + "_genomescope/summary.txt"
    seqstats = "RESULTS/STATISTICS/SEQKIT/" + PREFIX + ".ccs.readStats.txt"
    # assembly stats
    seqstatsassembl = expand("RESULTS/STATISTICS/SEQKIT/{basename}.assemblyStat.txt", basename = ASSEMBLY_BASENAMES)
    gfastats = expand("RESULTS/STATISTICS/GFASTATS/{basename}_gfastats_report.txt", basename = ASSEMBLY_BASENAMES)
    quast = expand("RESULTS/STATISTICS/QUAST/{basename}/quast.log", basename = ASSEMBLY_BASENAMES)
    gtseqstat = expand("RESULTS/STATISTICS/GTSEQSTATS/{basename}_gtseqstat.txt", basename = ASSEMBLY_BASENAMES)
    merqury1 = expand("RESULTS/STATISTICS/FASTK/{basename}_hifiasm.qv", basename=HIFIASM_MERQURY_BASENAMES)
    if DO_DEDUP: merqury2 = expand("RESULTS/STATISTICS/FASTK/{basename}_purge_dups.qv", basename = HIFIASM_MERQURY_BASENAMES)
    if DO_HIC: merqury3 = expand("RESULTS/STATISTICS/FASTK/{basename}_yahs.qv", basename = HIFIASM_MERQURY_BASENAMES)
    
    OUTPUT_FILES = ( 
    OUTPUT_FILES + 
    [katgc1] + [katgc2] + [katgc3] + 
    [genescope1] + 
    [seqstats] + seqstatsassembl + 
    quast + gfastats + gtseqstat + merqury1 + [smudgeplot]
    )
    if DO_DEDUP: OUTPUT_FILES = OUTPUT_FILES + merqury2
    if DO_HIC: OUTPUT_FILES = OUTPUT_FILES + merqury3

if DO_BUSCO:
    buscoout = expand("RESULTS/BUSCO/{basename}/run_{db}/missing_busco_list.tsv", db = config["busco"]["lineage"], basename = ASSEMBLY_BASENAMES)
    summary=expand("RESULTS/BUSCO/{basename}/short_summary.specific.{db}.{basename}.txt", db = config["busco"]["lineage"], basename = ASSEMBLY_BASENAMES)
    plot="RESULTS/BUSCO/busco_figure.png"
    OUTPUT_FILES = OUTPUT_FILES + buscoout + summary + [plot]
if DO_HIC:
    # yahs
    combined1 = expand("RESULTS/HIC/RAW/{basename}_{hic}.bam", hic = HIC_BASENAMES, basename = HIFIASM_BASENAMES)
    filt1 = expand("RESULTS/HIC/FILTERED/{basename}_{hic}_filtered.bam", hic = HIC_BASENAMES, basename = HIFIASM_BASENAMES)
    flagstat=expand("RESULTS/HIC/COMBINED/DEDUPED/{basename}_mapped_hic_reads_withRG_deduped_resorted_flagstat.tsv",basename = HIFIASM_BASENAMES)
    qualimap=expand("RESULTS/STATISTICS/QUALIMAP/{basename}/report.pdf",basename = HIFIASM_BASENAMES)
    #contact_map=expand("RESULTS/STATISTICS/HICSTUFF/{scaff_basename}/{scaff_basename}_map.png", scaff_basename = SCAF_BASENAMES)
    pretext=expand("RESULTS/STATISTICS/PRETEXT/{basename}.pretext",basename=HIFIASM_BASENAMES)
    OUTPUT_FILES  = OUTPUT_FILES + combined1 + filt1 + flagstat + qualimap + pretext # + contact_map
if DO_DEDUP:
    getseqs1 = expand("RESULTS/PURGE_DUPS/{basename}_seqs_purged.fasta", basename = HIFIASM_BASENAMES)
    getseqs2 = expand("RESULTS/MINIMAP2/{basename}_reads_split.paf", basename = HIFIASM_BASENAMES)
    getseqs3 = expand("RESULTS/MINIMAP2/{basename}_genome_split.paf", basename = HIFIASM_BASENAMES)
    getseqs4 = expand("RESULTS/PURGE_DUPS/{basename}_split.fasta", basename = HIFIASM_BASENAMES)

    OUTPUT_FILES =( OUTPUT_FILES + getseqs1 + getseqs2 + getseqs3 + getseqs4)
if DO_CCS:
    #checkfinished="RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/check.txt"
    # ccsreads="DATA/{basename}.subreads.bam"
    # hifireads="DATA/{basename}.fastq.gz"
    # ccspbi="DATA/{basename}.subreads.bam.pbi"
    ccsout = expand("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam", chunknumber = CHUNK_NMB)
    classification = "DATA/DECONTAMINATED/" + PREFIX + "_kraken_classification.txt"
    kraken_report = "DATA/DECONTAMINATED/" + PREFIX + "_kraken.report.csv"
    decontaminated = "DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz"
    contamination = "DATA/DECONTAMINATED/" + PREFIX + "_kraken_classified.fq.gz"
    OUTPUT_FILES = (OUTPUT_FILES + 
    ccsout  +  
    [classification] + [kraken_report] + [decontaminated] + [contamination] )
    # + [checkfinished] + [ccsreads] + [ccspbi] + [hifireads1]
elif not DO_CCS:
    # hifireads = HIFI if DO_ADAPTER_FILT and HIFI else "DATA/" + PREFIX + ".fastq.gz"
    classification = "DATA/DECONTAMINATED/" + PREFIX + "_kraken_classification.txt"
    kraken_report = "DATA/DECONTAMINATED/" + PREFIX + "_kraken.report.csv"
    decontaminated = "DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq.gz"
    contamination = "DATA/DECONTAMINATED/" + PREFIX + "_kraken_classified.fq.gz"
    OUTPUT_FILES = OUTPUT_FILES + [classification] + [kraken_report] + [decontaminated] + [contamination] # + hifireads

### rule all ###

rule all:
    input:
        OUTPUT_FILES


### Functions ### 

# Used for rules that need to run on files from different directories and with different basenames
# eg stats, scaffolding
def get_input(wildcards, files, basenames, wildcard_name):
    """
    This function returns the file path from a list of files by matching the file path to the current wildcard.
    The wildcard is estimated by snakemake based on the rule all.
    
    Parameters:
    wildcards (object): The wildcards object from Snakemake
    files (list): A list of file paths
    basenames (list): A list of basenames corresponding to the file paths.
    
    Example input:
    FILE_PATHS = [
    "path/to/file_1.txt",
    "path/to/file_2.txt",
    "path/to/file_n.txt"
    ]
    BASENAMES = [
        "file_1",
        "file_2",
        "file_n"
    ]
    Returns:
    str: The file path that matches the wildcard.
    Example use: 
    lambda wildcards: get_input(wildcards, FILE_PATHS, BASENAMES, 'my_wildcard_name')
    """
    # For each pair, it checks if the basename matches the wildcard. If it does, it yields the file path.
    # within its parentheses. If no item satisfies the condition, a StopIteration exception is raised.

    return next(
        (f for f, b in zip(files, basenames) if b == getattr(wildcards, wildcard_name))
    )

def get_input_hifiasm(wildcards, files, basenames, wildcard_name):
    basename = getattr(wildcards, wildcard_name)
    primary_suffixes = ["_primary.fasta"]
    alternate_suffixes = ["_alternate.fasta"]

    primary_file = next((f for f, b in zip(files, cycle(basenames)) if b == basename and any(f.endswith(suffix) for suffix in primary_suffixes)), None)
    alternate_file = next((f for f, b in zip(files, cycle(basenames)) if b == basename and any(f.endswith(suffix) for suffix in alternate_suffixes)), None)

    return primary_file, alternate_file

# used for gt seqstat and gfastats
def extract_max_genome_haploid_length(filename):
    """
    Extracts the maximum genome haploid length from a genescope summary file.

    This function reads the provided file line-by-line, searching for a line that contains 
    the string "Genome Haploid Length". Then, it extracts the maximum haploid length value 
    from that line and returns it. 
    If the value cannot be extracted or converted to an integer, it returns None.

    Parameters:
    - filename (str): Path to the genescope summary.
    Returns:
    - int: The maximum genome haploid length
    """
    try:
        with open(filename, 'r') as f:
            for line in f:
                if "Genome Haploid Length" in line:
                    lengths = [s.replace(',', '').strip() for s in line.split('bp') if s.strip()]
                    if len(lengths) > 1:
                        try:
                            return int(lengths[1])
                        except ValueError:
                            return None
        return None
    except:
        return None

def calculate_cov(max_length_file, genome_size_file):
    """
    Calculate approximate coverage of hifi reads for use in HyPo polishing

    Parameters:
    - max_length_file (str): Path to the file containing the maximum sequence length. The max length should be 
                             the 8th tab-separated value in the file (seqkit stats output)
    - genome_size_file (str): Path to the file containing the total genome size. The genome size should be 
                              the only integer value in the file (created by function above)

    Returns:
    - int: approximate coverage
    """
    try:
        with open(max_length_file, 'r') as f:
            max_length = int(f.readline().split('\t')[7])
        with open(genome_size_file, 'r') as f:
            genome_size = int(f.readline().strip())
        if max_length <= 0 or genome_size <= 0:
            raise ValueError("Max length and genome size must be positive integers.")
        return int(max_length / genome_size)
    except:
        print("404 error")



def hifiasm_options(basename):
    """
    Split the basename wildcard to get 
    the general and dup purge options for hifiasm.

    Parameters:
    - basename (object): Wildcard for the current run
    Returns:
    - options: String of all hifiasm options to use
    """
    keywords = basename.rsplit('_', 1)
    key, purge = keywords

    # Check if the key is in HIFIASM_OPTIONS
    if key in HIFIASM_OPTIONS:
        base_options = HIFIASM_OPTIONS[key]
    else:
        raise ValueError(f"Invalid key: {key}")

    # Format the purge options based on the name
    if purge == "l0":
        purge_options = "-l0"
    elif purge == "l1":
        purge_options = "-l1"
    elif purge == "l2":
        purge_options = "-l2"
    elif purge == "l3":
        purge_options = "-l3"
    else:
        raise ValueError(f"Invalid purge: {purge}")

    # Combine general options with purge options
    options = f"{base_options} {purge_options}"
    return options

def get_input_hifiasm(wildcards, files, basenames, wildcard_name):
    basename = getattr(wildcards, wildcard_name)  # This gets the current basename from wildcards
    primary_suffix = "_primary.fasta"
    alternate_suffix = "_alternate.fasta"

    # We filter the files before finding primary and alternate to simplify the search
    filtered_files = [f for f in files if basename in f]  # Filter files to only those that include the basename

    primary_file = next((f for f in filtered_files if f.endswith(primary_suffix)), None)
    alternate_file = next((f for f in filtered_files if f.endswith(alternate_suffix)), None)

    return primary_file, alternate_file


def get_input_purgedups(wildcards, files, basenames, wildcard_name):
    basename = getattr(wildcards, wildcard_name)
    primary_suffix = "_primary_seqs_purged.fasta"
    alternate_suffix = "_alternate_seqs_purged.fasta"

    # We filter the files before finding primary and alternate to simplify the search
    filtered_files = [f for f in files if basename in f]  # Filter files to only those that include the basename

    primary_file = next((f for f in filtered_files if f.endswith(primary_suffix)), None)
    alternate_file = next((f for f in filtered_files if f.endswith(alternate_suffix)), None)

    return primary_file, alternate_file

def get_input_scaff(wildcards, files, basenames, wildcard_name):
    basename = getattr(wildcards, wildcard_name)
    primary_suffix = "_primary_yahs_scaffolds_final.fa"
    alternate_suffix = "_alternate_yahs_scaffolds_final.fa"

    # We filter the files before finding primary and alternate to simplify the search
    filtered_files = [f for f in files if basename in f]  # Filter files to only those that include the basename

    primary_file = next((f for f in filtered_files if f.endswith(primary_suffix)), None)
    alternate_file = next((f for f in filtered_files if f.endswith(alternate_suffix)), None)

    return primary_file, alternate_file


include: "makefiles/CCS_DEEPCONSENSUS.smk"
include: "makefiles/STATS.smk"
include: "makefiles/MERQURY.smk"
include: "makefiles/DECONTAM.smk"
include: "makefiles/HIFIASM.smk"
include: "makefiles/BUSCO.smk"
include: "makefiles/PURGEDUPS.smk"
include: "makefiles/HIC.smk"

# Snakemake Pipeline for HiFi Genome Assemblies
# This pipeline processes HiFi Genome Assemblies

# Note:
# - Make sure the config.yaml file is properly configured

### Load Configurations ###

import re
import os
configfile: "config.yaml"

############  SET UP PATH TO INPUT FILES AND FILE NAME STEMS AND SAVE THEM AS WILDCARDS  ############
# subreads and fasta hifi reads must be a single file.
# hic reads need to be 2 files, R1 and R2

SUBREADS=config["SUBREADS"]
HIFI=config["HIFI"]
HIC1=config["HIC1"]
HIC2=config["HIC2"]
ENZYMES=config["ENZYMES"]
STOP_AT_HIFIASM=config["STOP_AT_HIFIASM"]
DO_STATS=config["DO_STATS"]
DO_BUSCO=config["DO_BUSCO"]
DO_HIC=config["DO_HIC"]
DO_CCS=config["DO_CCS"]
DO_DEDUP=config["DO_DEDUP"]
DO_ADAPT_FILT=config["DO_ADAPT_FILT"]
DO_CONTACT_MAP=config["DO_CONTACT_MAP"]
CHUNK_NMB=list(range(1,config["CHUNKS"]+1))

### Constants ###

# remove potential special characters and dots from prefix but keep underscores
# allow for completely unknown file names
PREFIX=re.sub('[^a-zA-Z0-9 \_]', '', config["FILE_PREFIX"])

# Used for the RG group in bwa when mapping the hic to the assembly
# subreads and fasta hifi reads must be a single file.
# TODO add a way to handle complex file names
# save anything after split('.')[0] as variable extension
if SUBREADS:
    #READ_PATH=os.path.dirname(SUBREADS)
    SAMPLE_LABEL=os.path.basename(SUBREADS).split('.')[0]
elif HIFI:
    #    print("Subreads not provided.")
    #READ_PATH=os.path.dirname(HIFI)
    SAMPLE_LABEL=os.path.basename(HIFI).split('.')[0]

############ Fasta and Read Files to Run Full Stats on:

# create lists of files depending on which step of the pipeline needs them.
# for example, purge_dups only needs the hifiasm output, but yahs needs the purge_dups output
# and assembly statistics like seqkit need all three assemblies (hifiasm, purgedups, yahs scaffolded)

## assemblies
ASSEMBLY_FILES=["RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_primary.fasta"] + ["RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_alternate.fasta"]

if DO_DEDUP:
    ASSEMBLY_FILES=ASSEMBLY_FILES + ["RESULTS/PURGE_DUPS/" + PREFIX + "_primary_seqs_purged.hap.fa"] + ["RESULTS/PURGE_DUPS/" + PREFIX + "_alternate_seqs_purged.hap.fa"]
if DO_HIC:
    ASSEMBLY_FILES=ASSEMBLY_FILES + ["RESULTS/HIC/YAHS/" + PREFIX + "_primary_yahs_scaffolds_final.fa"] + ["RESULTS/HIC/YAHS/" + PREFIX + "_alternate_yahs_scaffolds_final.fa"]

# Extracting basenames from ASSEMBLY_FILES
# these will be used to identify the different statistics used on each file
# later we define a function that matches those to the full path and uses the basename as a wildcard to name the output
ASSEMBLY_BASENAMES = [os.path.basename(f).split('.')[0] for f in ASSEMBLY_FILES]

# for rules that only need to run on the hifi assembl only
HIFIASM_FILES = ["RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_primary.fasta"] + ["RESULTS/GENOME_ASSEMBLY/" + PREFIX + "_alternate.fasta"]
HIFIASM_BASENAMES=[os.path.basename(f).split('.')[0] for f in HIFIASM_FILES]

SCAF_FILES = ["RESULTS/HIC/YAHS/" + PREFIX + "_primary_yahs_scaffolds_final.fa"]+["RESULTS/HIC/YAHS/" + PREFIX + "_alternate_yahs_scaffolds_final.fa"]
SCAF_BASENAMES=[os.path.basename(f).split('.')[0] for f in SCAF_FILES]

## Do the same for the read files

# reads will be decontaminated regardless of whether they are produced 
# from the subreads during the run or if 
# they are provided in the config file

HIFI_READ_FILES=["DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq"]

# need before- and after- file for seqkit stats
if HIFI: 
    HIFI_READ_FILES = HIFI_READ_FILES + [HIFI]

# if we give subreads as input, we still want to run stats on the hifi that get created during the run
if SUBREADS: 
    HIFI_READ_FILES = HIFI_READ_FILES + ["DATA/" + PREFIX + ".fastq.gz"]

HIC_READ_FILES=[]
if HIC1 and HIC2: 
    HIC_READ_FILES = HIC_READ_FILES + [HIC1] + [HIC2]

# Extracting basenames from HIFI_READ_FILES & HIC_READ_FILES
# these will be used to identify the different statistics used on each file
HIFI_BASENAMES = [os.path.basename(f).split('.')[0] for f in HIFI_READ_FILES]
HIC_BASENAMES = [os.path.basename(f).split('.')[0] for f in HIC_READ_FILES]


###  Rule All Files ###

# Add output files to rule all input, 
# so that the rules containing them will be included in the workflow.
# Create an empty list and append output file names to it depending on what the user wishes to do.
# If a pipeline rule is defined and its output isn't found on the all rule, it simply won't run

### Wildcard Constraints ###
# need to constrain wildcards so that smk can better resolve ambiguous wildcards and files
# this would not be necessary if we provided a regex for our input files but here we want to

hic_pattern = "|".join(HIC_BASENAMES)

wildcard_constraints:
    chunknumber="[0-9]+",
    hic=hic_pattern,
    basename=".*_.*"
    
### Set up input for rule all ###

OUTPUT_FILES = []
 
# the goal here is to conditionally add files to the list of outputs so that they match the operations we wish to run
if STOP_AT_HIFIASM:
    # make sure that if the user forgot to set all other routes to false that they don't run
    # then add the hifiasm output files to rule all, if we do other steps we wouldn't need these files in there
    DO_STATS = False
    DO_BUSCO = False
    DO_HIC = False
    DO_DEDUP = False
    DO_ADAPT_FILT = False
    hifiasm = expand("RESULTS/GENOME_ASSEMBLY/{basename}.fasta", basename=HIFIASM_BASENAMES)

    OUTPUT_FILES=OUTPUT_FILES + hifiasm
#assembly stats for all assemblies
if DO_STATS:

    # read stats
    katgc1="RESULTS/STATISTICS/FASTK/" + PREFIX + ".st.png"
    katgc2="RESULTS/STATISTICS/FASTK/" + PREFIX + ".fi.png"
    katgc3="RESULTS/STATISTICS/FASTK/" + PREFIX + ".ln.png"
    ploidyplot1="RESULTS/STATISTICS/FASTK/PLOIDYPLOT/" + PREFIX + ".png"
    genescope1="RESULTS/STATISTICS/FASTK/genescopeFK/" + PREFIX + "_progress.txt"
    seqstats="RESULTS/STATISTICS/SEQKIT/" + PREFIX + ".ccs.readStats.txt"

    # assembly stats
    seqstatsassembl=expand("RESULTS/STATISTICS/SEQKIT/{basename}.assemblyStat.txt", basename=ASSEMBLY_BASENAMES)
    gfastats=expand("RESULTS/STATISTICS/GFASTATS/{basename}_gfastats_report.txt", basename=ASSEMBLY_BASENAMES)
    gaasstats=expand("RESULTS/STATISTICS/GAAS/{basename}_gaasstats.txt", basename=ASSEMBLY_BASENAMES)
    quast=expand("RESULTS/STATISTICS/{basename}_hifiasm/quast/quast.log", basename=ASSEMBLY_BASENAMES)
    gtseqstat=expand("RESULTS/STATISTICS/GTSEQSTATS/{basename}_gtseqstat.txt", basename=ASSEMBLY_BASENAMES)
    merqury=expand("RESULTS/STATISTICS/FASTK/{basename}.completeness.stats", basename=ASSEMBLY_BASENAMES)


    OUTPUT_FILES = ( 
    OUTPUT_FILES + 
    [katgc1] + [katgc2] + [katgc3] + 
    [ploidyplot1] + [genescope1] + 
    [seqstats] + seqstatsassembl + 
    quast + gfastats + gaasstats + gtseqstat + [merqury]
    )

#regular hifiasm busco for both haplos
if DO_BUSCO:
    buscoout=expand("RESULTS/BUSCO/{basename}/{basename}/run_{db}/missing_busco_list.tsv", db=config["LINEAGE"], basename=HIFIASM_BASENAMES)
    OUTPUT_FILES=OUTPUT_FILES + buscoout

# turn path string to list and add to output
if DO_HIC:
    hicresorted1=expand("RESULTS/HIC/YAHS/{basename}.fa", basename=SCAF_BASENAMES)
    combined1=expand("RESULTS/HIC/RAW/{basename}_{hic}.bam", hic=HIC_BASENAMES, basename=HIFIASM_BASENAMES)
    filt1=expand("RESULTS/HIC/FILTERED/{basename}_{hic}_filtered.bam", hic=HIC_BASENAMES, basename=HIFIASM_BASENAMES)
    busco_scaff=expand("RESULTS/BUSCO/{basename}/{basename}/run_{db}/missing_busco_list.tsv", db=config["LINEAGE"], basename=SCAF_BASENAMES)
    OUTPUT_FILES=OUTPUT_FILES + [hicresorted1] + combined1 + filt1 + [busco_scaff]

# do purge_dups for both haplotypes
if DO_DEDUP:
    buscodedup=expand("RESULTS/BUSCO/{basename}_seqs_purged/{basename}_seqs_purged/run_{db}/missing_busco_list.tsv", db=config["LINEAGE"], basename=HIFIASM_BASENAMES)
    getseqs1=expand("RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa", basename=HIFIASM_BASENAMES)
    getseqs2=expand("RESULTS/MINIMAP2/{basename}_reads_split.paf", basename=HIFIASM_BASENAMES)
    getseqs3=expand("RESULTS/MINIMAP2/{basename}_genome_split.paf", basename=HIFIASM_BASENAMES)
    getseqs4=expand("RESULTS/PURGE_DUPS/{basename}_split.fasta", basename=HIFIASM_BASENAMES)
    OUTPUT_FILES=OUTPUT_FILES + getseqs1 + buscodedup + getseqs2 + getseqs3 + getseqs4

# perform decontamination regardless of whether the user provides bam subreads or hifi fastq reads
if DO_CCS:
    #checkfinished="RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/check.txt"
    # ccsreads="DATA/{basename}.subreads.bam"
    # hifireads="DATA/{basename}.fastq.gz"
    # ccspbi="DATA/{basename}.subreads.bam.pbi"
    mergefq="RESULTS/PREPROCESSING/CCS_PACBIO/MERGED/" + PREFIX + ".merged.ccs.fastq"
    ccsout=expand("RESULTS/PREPROCESSING/CCS_PACBIO/" + PREFIX + "_{chunknumber}.ccs.bam", chunknumber=CHUNK_NMB)
    classification="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classification.txt"
    kraken_report="DATA/DECONTAMINATED/" + PREFIX + "_kraken.report.csv"
    decontaminated="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq"
    contamination="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classified.fq"

    OUTPUT_FILES = (OUTPUT_FILES + 
    ccsout  + [mergefq]  + 
    [classification] + [kraken_report] + [decontaminated] + [contamination] )
    # + [checkfinished] + [ccsreads] + [ccspbi] + [hifireads1]

elif not DO_CCS:
    # hifireads = HIFI if DO_ADAPTER_FILT and HIFI else "DATA/" + PREFIX + ".fastq.gz"
    classification="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classification.txt"
    kraken_report="DATA/DECONTAMINATED/" + PREFIX + "_kraken.report.csv"
    decontaminated="DATA/DECONTAMINATED/" + PREFIX + "_kraken_unclassified.fq"
    contamination="DATA/DECONTAMINATED/" + PREFIX + "_kraken_classified.fq"
    OUTPUT_FILES=OUTPUT_FILES + [classification] + [kraken_report] + [decontaminated] + [contamination] # + hifireads1
if DO_CONTACT_MAP:
    nobinning=expand("RESULTS/HIC_CONTACT_MAP_{scaff_basename}/{scaff_basename}_map.png", scaff_basename=SCAF_BASENAMES)
    fivekbinning=expand("RESULTS/HIC_CONTACT_MAP_{scaff_basename}/{scaff_basename}_5kmap.png", scaff_basename=SCAF_BASENAMES)
    tenkbinning=expand("RESULTS/HIC_CONTACT_MAP_{scaff_basename}/{scaff_basename}_10kmap.png", scaff_basename=SCAF_BASENAMES)
    OUTPUT_FILES=OUTPUT_FILES + nobinning + fivekbinning + tenkbinning

### Rule All ###

# Terminal node:
rule all:
    input:
        OUTPUT_FILES


### Functions ### 

# Used for rules that need to run on files from different directories and with different basenames
# eg stats, scaffolding
def get_input(wildcards, files, basenames, wildcard_name):
    """
    This function returns the file path from a list of files based on a given wildcard.
    
    Parameters:
    wildcards (object): The wildcards object from Snakemake, which contains the wildcard values.
    files (list): A list of file paths.
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

# used for gt seqstat and gfastats
def extract_max_genome_haploid_length(filename):
    """
    Extracts the maximum genome haploid length from a genescope summary file.

    This function reads the provided file line-by-line, searching for a line that contains 
    the string "Genome Haploid Length". Then, it extracts the maximum haploid length value 
    from that line and returns it. 
    If the value cannot be extracted or converted to an integer, it returns None.

    Parameters:
    - filename (str): Relative or absolute path to the genescope summary.

    Returns:
    - int: The maximum genome haploid length if found and successfully extracted.
    - None: If the length cannot be found, extracted, or converted to an integer.
    """
    with open(filename, 'r') as f:
        for line in f:
            if "Genome Haploid Length" in line:
                # Extract the max value safely
                lengths = [s.replace(',', '').strip() for s in line.split('bp') if s.strip()]
                if len(lengths) > 1:
                    try:
                        return int(lengths[1])
                    except ValueError:
                        return None
    return None


include: "makefiles/CCS_DEEPCONSENSUS.smk"
include: "makefiles/STATS.smk"
include: "makefiles/DECONTAM.smk"
include: "makefiles/HIFIASM.smk"
include: "makefiles/BUSCO.smk"
include: "makefiles/PURGEDUPS.smk"
include: "makefiles/HIC.smk"
include: "makefiles/HICSTUFF_COOL.smk"

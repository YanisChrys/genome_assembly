# Step 1: Chromosome Extraction
# Extracts individual chromosomes using `awk`.

# Step 2: Chromosome Indexing
# Indexes each chromosome file with `samtools`.

# Step 3: Bed File Creation
# Creates a BED file for each chromosome using `awk`.

# Step 4: Reads Retrieval
# Retrieves the reads for each SMRT cell in FASTQ format.

# Step 5: Reads Alignment
# Aligns reads to the reference genome using `pbmm2`, outputting a BAM file.

# Step 6: Alignment Splitting
# Splits the aligned reads by chromosome with `samtools`, creating a separate BAM file per chromosome.

# Step 7: BAM Files Merging
# Merges the BAM files for each SMRT cell into a single BAM file per chromosome using `samtools`.

# Step 8: Variant Calling
# Calls variants from the aligned reads using `DeepVariant`.

# Step 9: Variant Filtering
# Filters variants using `bcftools`, retaining only those that pass quality filters.

# Step 10: Consensus Sequence Creation
# Creates a consensus sequence for each chromosome using `bcftools`.

# Step 11: Quality Assessment
# Assesses the quality of the consensus sequence using `Merqury`.


"""
Snakemake Pipeline for HiFi Genome Assemblies
This pipeline processes HiFi Genome Assemblies
             
Version: 1.0
Author: Ioannis Chrysostomakis
Email: john.chrysost@gmail.com
Date Created: YYYY-MM-DD
Last Modified: YYYY-MM-DD

License: MIT

Usage: 

Dependencies: 

Acknowledgements: 

Notes: Make sure the configuration file is properly set up 
and make sure to read the README file before attempting to run this pipeline
"""

import re
import os
import glob
configfile: "config/config_polish_deepvariant.yaml"

# get sample basename regardless of the extension (fa, fasta) and create sample wildcard
asm_path=config["asm_path"]
extensions = ["fa", "fasta","fa.gz", "fasta.gz"]
patterns = [os.path.join(asm_path, f"*.{ext}") for ext in extensions] # create list of possible targets
files = [file for pattern in patterns for file in glob.glob(pattern)] # now look up these combinations and find the true files
samples = {os.path.basename(file).rsplit('.', 1)[0] for file in files} # create the sample wildcard out of those true samples

# get the basename of the reads so fastk only runs once
PREFIX = re.sub(r'\.(fq\.gz|fastq\.gz|fq|fastq)$', '', os.path.basename(config["reads"]))

# characterist of your run to be added to outputs:
run=config["run"]

# wildcard_constraints:
#     basename=".*_.*",
#     chrom = "[a-zA-Z0-9_-]*"

asm_path=config["asm_path"]
extensions = ["fa", "fasta","fa.gz", "fasta.gz"]
patterns = [os.path.join(asm_path, f"*.{ext}") for ext in extensions] # create list of possible targets
files = [file for pattern in patterns for file in glob.glob(pattern)] # now look up these combinations and find the true files
samples = {os.path.basename(file).rsplit('.', 1)[0] for file in files}

rule all: 
    input:
        expand("RESULTS/ASM_POLISHING_DEEPVARIANT/{sample}_" + run + "/interval.list",sample=samples)



# Function:

# def get_chromosomes():
# 	with open('ref_genom/interval.list', 'r') as file:
# 		return file.read().splitlines()


# chromosomes = get_chromosomes()


#replace commas and spaces with underscores
rule correct_headers:
    input:
        lambda wildcards: next(os.path.join(asm_path, f"{wildcards.sample}.{ext}")
                               for ext in extensions
                               if glob.glob(os.path.join(asm_path, f"{wildcards.sample}.{ext}")))
    output:
        "RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta"
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit replace -p '[ ,]' -r '_' {input} > {output}
    """

rule index_ref:
    input:
        "RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta"
    output:
        fai="RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta.fai"
    conda:
        "../envs/samtools_minimap.yaml"
    shell: """
        samtools faidx {input}
    """


# split genome into chromosomes
# two-pass: read file twice to reduce memory usage
# split by the full id name of each chromosome
checkpoint split_by_chromosomes:
    input:
        fasta="RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta",
        fai="RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta.fai"
    output:
        directory("RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_" + run + "/")
    conda:
        "../envs/seqkit.yaml"
    params:
        outdir="RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_" + run + "/"
    shell: """
        seqkit split --quiet --two-pass -i --by-id-prefix "" --out-dir {params.outdir} {input.fasta}
    """

# delete if not useful
rule get_chromosomes:
    input: 
        fai="RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta.fai",
        dir="RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_" + run + "/"
    output:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/{sample}_" + run + "/interval.list"
    params:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/{sample}_" + run + "/"
    shell: """
        awk '{{print $1}}' {input.fai} > {output}
    """


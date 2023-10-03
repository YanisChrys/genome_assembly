import os
import re

# path to hic without final backslash
PATH_TO_RAW_HIC="/share/pool/CompGenomVert/RawData/UrPie/DATA"
PATH_TO_ADAPTERS='/share/pool/CompGenomVert/RawData/UrPie/adapters.fasta'
TEMP_DIR="/share/pool/ychrysostomakis/tmp"
KRAKEN_DB="/share/pool/databases/kraken"

SAMPLES = glob_wildcards(PATH_TO_RAW_HIC + "/{sample}_R1.fastq.gz").sample

rule all:
    input:
    #FASTQC MULTIQC
        expand("output/fastqc/{sample}_corrected_R1_fastqc.zip", sample=SAMPLES),
        expand("output/fastqc/{sample}_corrected_R2_fastqc.zip", sample=SAMPLES),
        "output/multiqc/multiqc_report.html"

### Functions ###

def get_input_path(sample, type):
    """
    Returns file paths based on the raw input basename and a unique command identifier.

    Parameters:
    - sample (str): The sample name.
    - type (str): The type of file (raw, trimmed, corrected, decontaminated).

    Returns:
    - tuple: A tuple containing paths for R1 and R2.
    """

    if type == "raw":
        return (
            f"{PATH_TO_RAW_HIC}/{sample}_R1.fastq.gz",
            f"{PATH_TO_RAW_HIC}/{sample}_R2.fastq.gz"
        )
    elif type == "trimmed":
        return (
            f"output/qc/{sample}_trimmed_R1.fastq.gz",
            f"output/qc/{sample}_trimmed_R2.fastq.gz"
        )
    elif type == "corrected":
        return (
            f"output/{sample}_corrected_R1.fastq.gz",
            f"output/{sample}_corrected_R2.fastq.gz"
        )
    elif type == "decontaminated":
        return (
            f"output/kraken2/{sample}_decontaminated_R1.fastq",
            f"output/kraken2/{sample}_decontaminated_R2.fastq"
        )
    else:
        raise ValueError(f"Unknown type: {type}")

### Rules ###

# Trimmomatic Adapter and Quality Trimming:
# - Assumes input quality scores are in Phred+33 format.
# - Generates a log file detailing the trimming process.
# - 'ILLUMINACLIP' removes Illumina adapters. Parameters are:
#   - Adapter sequence path.
#   - Seed mismatches (2): specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
#   - Simple clip threshold (30): specifies how accurate the match between any adapter etc. sequence must be against a read.
#   - Min adapter length (10): In addition to the alignment score, palindrome mode can verify that a minimum length of adapter has been detected.
# - 'LEADING' and 'TRAILING' remove bases from the start or end of a read
# - 'SLIDINGWINDOW' starts scanning at the 5' end and clips the read once the average quality within the window falls below the threshold
# - 'MINLEN' drops the read if it is below a specified length


rule TRIMMOMATIC:
    input:
        r1= PATH_TO_RAW_HIC + '/{sample}_R1.fastq.gz',
        r2= PATH_TO_RAW_HIC + '/{sample}_R2.fastq.gz'
    output:
        r1_paired='output/qc/{sample}_trimmed_R1.fastq.gz',
        r2_paired='output/qc/{sample}_trimmed_R2.fastq.gz',
        r1_unpaired="output/qc/{sample}_unpaired_R1.fastq.gz",
        r2_unpaired="output/qc/{sample}_unpaired_R2.fastq.gz",
        log="output/qc/trimlog_{sample}.log"
    params:
        adapter_path=PATH_TO_ADAPTERS
    threads:
        workflow.cores
    shell: """
        module load trimmomatic/0.39

        java -jar $TRIMMOMATIC/trimmomatic-0.39.jar PE \
        -threads {threads} -phred33 \
        -trimlog {output.log} \
        {input.r1} {input.r2} \
        {output.r1_paired} \
        {output.r1_unpaired} \
        {output.r2_paired} \
        {output.r2_unpaired} \
        ILLUMINACLIP:{params.adapter_path}:2:30:10 \
        LEADING:5 TRAILING:5 \
        SLIDINGWINDOW:5:20 MINLEN:30
    """

rule KRAKEN2:
    input:
        r1='output/qc/{sample}_trimmed_R1.fastq.gz',
        r2='output/qc/{sample}_trimmed_R2.fastq.gz'
    output:
        classification="output/kraken2/{sample}_kraken_classification.txt",
        report="output/kraken2/{sample}_kraken.report.csv",
        decontaminated1="output/kraken2/{sample}_decontaminated_R1.fastq",
        contamination1="output/kraken2/{sample}_kraken_classified_R1.fastq",
        decontaminated2="output/kraken2/{sample}_decontaminated_R2.fastq",
        contamination2="output/kraken2/{sample}_kraken_classified_R2.fastq"
    threads:
        workflow.cores
    params:
        unclassified_input="output/kraken2/{sample}_decontaminated_R#.fastq",
        classified_input="output/kraken2/{sample}_kraken_classified_R#.fastq",
        krakendb=KRAKEN_DB
    conda:
        "envs/kraken2.yaml"

    shell: """
        kraken2 --db {params.krakendb}
            --threads {threads} \
            --use-names \
            --confidence 0.51 \
            --output {output.classification} \
            --report {output.report} \
            --unclassified-out {params.unclassified_input} \
            --classified-out {params.classified_input} \
            {input.r1} {input.r2}
    """


# Tadpole Error Correction:
# - Uses multiple threads for parallel processing
# - Takes paired-end input reads (in1 and in2) for error correction
# - Outputs corrected reads (out1 and out2)
# - 'reassemble': enabling reassembly-based correction
# - 'mode' set to 'correct', so we do error correction
# - 'minprob': probability for making a correction
# - 'prefilter' and 'prehashes': speed up processing
# - 'prealloc' preallocates memory based on 'prefiltersize' to speed up processing
# - 'fillfast' accelerates the process.
# - 'pincer', 'ecc', and 'tail' are error-correction techniques enabled for the process
# - 'k': kmer length for error correction.
# - 'prefiltersize': fraction of memory preallocated for the prefilter

rule TADPOLE_ERROR_CORRECTION:
    input:
        r1='output/kraken2/{sample}_decontaminated_R1.fastq',
        r2='output/kraken2/{sample}_decontaminated_R2.fastq'
    output:
        r1_corrected='output/{sample}_corrected_R1.fastq.gz',
        r2_corrected='output/{sample}_corrected_R2.fastq.gz'
    threads:
        workflow.cores
    conda:
        "envs/tadpole.yaml"
    shell: """
        tadpole.sh \
        threads={threads} \
        in1={input.r1} \
        in2={input.r2} \
        out1={output.r1_corrected} \
        out2={output.r2_corrected} \
        reassemble=t \
        mode=correct \
        minprob=0.6 \
        prefilter=1 \
        prehashes=2 \
        prealloc=t \
        fillfast=t \
        pincer=t \
        ecc=t \
        tail=t \
        k=50 \
        prefiltersize=0.2
        """

rule FASTQC:
    input:
        lambda wildcards: get_input_path(wildcards.sample, wildcards.type)
    output:
        r1_html='output/fastqc/{sample}_{type}_R1_fastqc.zip',
        r2_html='output/fastqc/{sample}_{type}_R2_fastqc.zip'
    params:
        outdir="output/fastqc/",
        tmpdir=TEMP_DIR
    threads:
        min(workflow.cores,6)
    conda:
        "envs/multiqc.yaml"
    priority: 1
    shell: """
        fastqc -o {params.outdir} -t {threads} -d {params.tmpdir} {input[0]} {input[1]}
    """


rule MULTIQC:
    input:
        expand('output/fastqc/{sample}_raw_R1_fastqc.zip', sample=SAMPLES),
        expand('output/fastqc/{sample}_raw_R2_fastqc.zip', sample=SAMPLES),
        expand('output/fastqc/{sample}_trimmed_R1_fastqc.zip', sample=SAMPLES),
        expand('output/fastqc/{sample}_trimmed_R2_fastqc.zip', sample=SAMPLES),
        expand('output/fastqc/{sample}_decontaminated_R1_fastqc.zip', sample=SAMPLES),
        expand('output/fastqc/{sample}_decontaminated_R1_fastqc.zip', sample=SAMPLES),
        expand('output/fastqc/{sample}_corrected_R1_fastqc.zip', sample=SAMPLES),
        expand('output/fastqc/{sample}_corrected_R2_fastqc.zip', sample=SAMPLES)
    output:
        "output/multiqc/multiqc_report.html"
    params:
        qcdir="output/fastqc/",
        outdir="output/multiqc/"
    conda:
        "envs/multiqc.yaml"
    shell: """
        multiqc {params.qcdir} -o {params.outdir} -p
    """
import os
import re
configfile: "config/config.yaml"
# path to hic without final backslash(important in rule paths)
PATH_TO_ADAPTERS="/home/ychrysostomakis/adapters/uniq_adapters.fasta"
TEMP_DIR="/share/pool/ychrysostomakis/tmp"
KRAKEN_DB="/share/pool/databases/kraken"
MIN_READ_LENGTH = 95
HIC1 = config["hic1"]
HIC2 = config["hic2"]
# return run prefix regardless of whether it ends in .R1 or _R1
PREFIX = re.search(r'(.+?)[._]R1', os.path.basename(HIC1)).group(1)





rule all:
    input:
    #FASTQC MULTIQC
        "RESULTS/HICREADS/MULTIQC/multiqc_report.html"

### Functions ###

# for when you have a wildcard
# def get_input_path(sample, type):
#     """
#     Returns file paths based on the raw input basename and a unique command identifier.

#     Parameters:
#     - sample (str): The sample name.
#     - type (str): The type of file (raw, trimmed, corrected, decontaminated).

#     Returns:
#     - tuple: A tuple containing paths for R1 and R2.
#     """

#     if type == "raw":
#         return (
#             f"{PATH_TO_RAW_HIC}/{sample}_R1.fastq.gz",
#             f"{PATH_TO_RAW_HIC}/{sample}_R2.fastq.gz"
#         )
#     elif type == "trimmed":
#         return (
#             f"output/qc/{sample}_trimmed_R1.fastq.gz",
#             f"output/qc/{sample}_trimmed_R2.fastq.gz"
#         )
#     elif type == "corrected":
#         return (
#             f"output/{sample}_corrected_R1.fastq.gz",
#             f"output/{sample}_corrected_R2.fastq.gz"
#         )
#     elif type == "decontaminated":
#         return (
#             f"output/kraken2/{sample}_decontaminated_R_1.fastq",
#             f"output/kraken2/{sample}_decontaminated_R_2.fastq"
#         )
#     else:
#         raise ValueError(f"Unknown type: {type}")

def get_input_path(type):
    """
    Returns file paths based on the PREFIX and a unique command identifier.

    Parameters:
    - type (str): The type of file (raw, trimmed, corrected, decontaminated).

    Returns:
    - tuple: A tuple containing paths for R1 and R2.
    """
    type_paths = {
        "_": (HIC1, HIC2),
        "_trimmed_": (f"RESULTS/HICREADS/QC/{PREFIX}_trimmed_R1.fastq.gz",
                    f"RESULTS/HICREADS/QC/{PREFIX}_trimmed_R2.fastq.gz"),
        "_corrected_": (f"RESULTS/HICREADS/QC/{PREFIX}_corrected_R1.fastq.gz",
                      f"RESULTS/HICREADS/QC/{PREFIX}_corrected_R2.fastq.gz")
    }

    if type in type_paths:
        return type_paths[type]
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


rule FASTP:
    input:
        r1=HIC1,
        r2=HIC2
    output:
        r1_paired="RESULTS/HICREADS/QC/" + PREFIX + "_trimmed_R1.fastq.gz",
        r2_paired="RESULTS/HICREADS/QC/" + PREFIX + "_trimmed_R2.fastq.gz",
        report_html="RESULTS/HICREADS/QC/fastp_30_20_report_" + PREFIX + ".html",
        report_json="RESULTS/HICREADS/QC/fastp_30_20_report_" + PREFIX + ".json"
    params:
        adapter_path=PATH_TO_ADAPTERS,
        min_len=MIN_READ_LENGTH
    threads:
        min(workflow.cores,10)
    conda:
        "../envs/fastp.yaml"
    shell: """
        fastp -V -w {threads} \
        --adapter_fasta {params.adapter_path} \
        --length_required {params.min_len} \
        --qualified_quality_phred 20 \
        --trim_front1 5 \
        --trim_front2 5 \
        -i {input.r1} -I {input.r2} \
        -o {output.r1_paired} -O {output.r2_paired} \
        -h {output.report_html} \
        -j {output.report_json}
    """

rule KRAKEN2:
    input:
        r1="RESULTS/HICREADS/QC/" + PREFIX + "_trimmed_R1.fastq.gz",
        r2="RESULTS/HICREADS/QC/" + PREFIX + "_trimmed_R2.fastq.gz"
    output:
        classification="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_kraken_classification.txt",
        report="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_kraken.report.csv",
        decontaminated1="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_decontaminated_R_1.fastq",
        contamination1="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_kraken_classified_R_1.fastq",
        decontaminated2="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_decontaminated_R_2.fastq",
        contamination2="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_kraken_classified_R_2.fastq"
    threads:
        min(workflow.cores,10)
    params:
        unclassified_input="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_decontaminated_R#.fastq",
        classified_input="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_kraken_classified_R#.fastq",
        krakendb=KRAKEN_DB
    conda:
        "../envs/kraken2.yaml"
    shell: """
        kraken2 --db {params.krakendb} \
            --threads {threads} \
            --use-names \
            --paired \
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
        r1="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_decontaminated_R_1.fastq",
        r2="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_decontaminated_R_2.fastq"
    output:
        r1_corrected="RESULTS/HICREADS/QC/" + PREFIX + "_corrected_R1.fastq.gz",
        r2_corrected="RESULTS/HICREADS/QC/" + PREFIX + "_corrected_R2.fastq.gz"
    threads:
        min(workflow.cores,10)
    conda:
        "../envs/tadpole.yaml"
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

rule FASTQC_DECONTAM:
    input:
        decontaminated1="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_decontaminated_R_1.fastq",
        decontaminated2="RESULTS/HICREADS/KRAKEN2/" + PREFIX + "_decontaminated_R_2.fastq",
    output:
        r1_html="RESULTS/HICREADS/FASTQC/" + PREFIX + "_decontaminated_R_1_fastqc.zip",
        r2_html="RESULTS/HICREADS/FASTQC/" + PREFIX + "_decontaminated_R_2_fastqc.zip"
    params:
        outdir="RESULTS/HICREADS/FASTQC/",
        tmpdir=TEMP_DIR
    threads:
        min(workflow.cores,10)
    conda:
        "../envs/multiqc.yaml"
    priority: 1
    shell: """
        fastqc -o {params.outdir} -t {threads} -d {params.tmpdir} {input.decontaminated1} {input.decontaminated2}
    """

rule FASTQC:
    input:
        lambda wildcards: get_input_path(wildcards.type)
    output:
        r1_html="RESULTS/HICREADS/FASTQC/" + PREFIX + "{type}R1_fastqc.zip",
        r2_html="RESULTS/HICREADS/FASTQC/" + PREFIX + "{type}R2_fastqc.zip"
    params:
        outdir="RESULTS/HICREADS/FASTQC/",
        tmpdir=TEMP_DIR
    threads:
        min(workflow.cores,10)
    conda:
        "../envs/multiqc.yaml"
    priority: 1
    shell: """
        fastqc -o {params.outdir} -t {threads} -d {params.tmpdir} {input[0]} {input[1]}
    """

rule MULTIQC:
    input:
        "RESULTS/HICREADS/FASTQC/" + PREFIX + "_R1_fastqc.zip",
        "RESULTS/HICREADS/FASTQC/" + PREFIX + "_R2_fastqc.zip",
        "RESULTS/HICREADS/FASTQC/" + PREFIX + "_trimmed_R1_fastqc.zip",
        "RESULTS/HICREADS/FASTQC/" + PREFIX + "_trimmed_R2_fastqc.zip",
        "RESULTS/HICREADS/FASTQC/" + PREFIX + "_decontaminated_R_1_fastqc.zip",
        "RESULTS/HICREADS/FASTQC/" + PREFIX + "_decontaminated_R_2_fastqc.zip",
        "RESULTS/HICREADS/FASTQC/" + PREFIX + "_corrected_R1_fastqc.zip",
        "RESULTS/HICREADS/FASTQC/" + PREFIX + "_corrected_R2_fastqc.zip"
    output:
        "RESULTS/HICREADS/MULTIQC/multiqc_report.html"
    params:
        qcdir="RESULTS/HICREADS/FASTQC/",
        outdir="RESULTS/HICREADS/MULTIQC/"
    conda:
        "../envs/multiqc.yaml"
    shell: """
        multiqc {params.qcdir} -o {params.outdir} -p
    """


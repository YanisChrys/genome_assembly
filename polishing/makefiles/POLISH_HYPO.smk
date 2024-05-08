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
configfile: "config/config_polish_hypo.yaml"

# get sample basename regardless of the extension (fa, fasta) and create sample wildcard
asm_path=config["asm_path"]

asm_path="/home/ychrysostomakis/EpiAus/asm_fasta/"
extensions = ["fa", "fasta","fa.gz", "fasta.gz"]
patterns = [os.path.join(asm_path, f"*.{ext}") for ext in extensions] # create list of possible targets
files = [file for pattern in patterns for file in glob.glob(pattern)] # now look up these combinations and find the true files
samples = {os.path.basename(file).rsplit('.', 1)[0] for file in files} # create the sample wildcard out of those true samples

samples=["EpiAus_202401_3SMRT_hifiasm_incpert_D10N150_l2_primary"]
# get the basename of the reads so fastk only runs once
PREFIX = re.sub(r'\.(fq\.gz|fastq\.gz|fq|fastq)$', '', os.path.basename(config["reads"]))

# characterist of your run to be added to outputs:
run=config["run"]

def get_input(wildcards, files, basenames, wildcard_name):
    return next(
        (f for f, b in zip(files, basenames) if b == getattr(wildcards, wildcard_name))
    )

FA_FILES = (expand("RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-1.fa", sample = samples) +
    expand("RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-2.fa", sample = samples) +
    expand("RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-3.fa", sample = samples))

BASENAMES = [os.path.basename(f).split('.')[0] for f in FA_FILES]


# wildcard_constraints:
#     basename=".*_.*",
#     chrom = "[a-zA-Z0-9_-]*"

rule all:
    input:
    # polishing
        expand("RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-3.fa", sample = samples),
    # stats
        expand("RESULTS/BUSCO/{basename}/run_{db}/missing_busco_list.tsv", db = config["busco"]["lineage"], basename=BASENAMES),
        expand("RESULTS/BUSCO/{basename}/short_summary.specific.{db}.{basename}.txt", db = config["busco"]["lineage"], basename=BASENAMES),
        expand("RESULTS/STATISTICS_POL/GFASTATS/{basename}_gfastats_report.txt", basename=BASENAMES),
        expand("RESULTS/STATISTICS_POL/QUAST/{basename}/quast.log", basename=BASENAMES),
        expand("RESULTS/STATISTICS_POL/FASTK/{basename}.completeness.stats", basename=BASENAMES),
        expand("RESULTS/STATISTICS_POL/SEQKIT/{basename}.assemblyStat.txt", basename=BASENAMES)

def calculate_cov(total_read_length, genome_size):
    """
    Calculate approximate coverage of hifi reads for use in HyPo polishing
    """
    total_read_length = float(total_read_length)
    genome_size = float(genome_size)
    
    # Perform division, truncate the result, and return coverage as an integer
    coverage = int(total_read_length/genome_size)
    return coverage

#replace commas and spaces with underscores
rule correct_headers:
    input:
        lambda wildcards: next(os.path.join(asm_path, f"{wildcards.sample}.{ext}")
                               for ext in extensions
                               if glob.glob(os.path.join(asm_path, f"{wildcards.sample}.{ext}")))
    output:
        touch("RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta")
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

rule run_hypo:
    input:
        asm="RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta",
        fai="RESULTS/CORRECTED_HEADERS/{sample}_corrected_headers.fasta.fai"
    output:
        out1="RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-1.fa",
        out2="RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-2.fa",
        out3="RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-3.fa",
        bam1=temp("RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-1.bam"),
        bam2=temp("RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-2.bam"),
        bam3=temp("RESULTS/ASM_POLISHING_HYPO/{sample}_" + run + "-3.bam")
    threads: 
        workflow.cores
    priority: 1
    params:
        genome_size=config["genome_size"],
        coverage = lambda wildcards: calculate_cov(config["total_read_length"],config["genome_size"]),
        pacbioreads=config["reads"],
        fofn="RESULTS/ASM_POLISHING_HYPO/{sample}_inputs.fofn"
    conda:
        "../envs/hypo.yaml"
    shell: """
        cat <<EOF > {params.fofn} 
{params.pacbioreads}
EOF

        minimap2 --secondary=no --MD --eqx -ax map-hifi -t {threads} {input.asm} $(echo {params.pacbioreads}) | samtools view -Sb - | samtools sort -o {output.bam1} -@ {threads} -
        samtools index {output.bam1}
        hypo -d {input.asm} -r @{params.fofn} -s {params.genome_size} -c {params.coverage} -b {output.bam1} -p 80 -t {threads} -o {output.out1}

        minimap2 --secondary=no --MD --eqx -ax map-hifi -t {threads} {output.out1} $(echo {params.pacbioreads}) | samtools view -Sb - | samtools sort -o {output.bam2} -@ {threads} -
        samtools index {output.bam2}
        hypo -d {output.out1} -r @{params.fofn} -s {params.genome_size} -c {params.coverage} -b {output.bam2} -p 80 -t {threads} -o {output.out2}

        minimap2 --secondary=no --MD --eqx -ax map-hifi -t {threads} {output.out2} $(echo {params.pacbioreads}) | samtools view -Sb - | samtools sort -o {output.bam3} -@ {threads} -
        samtools index {output.bam3}
        hypo -d {output.out2} -r @{params.fofn} -s {params.genome_size} -c {params.coverage} -b {output.bam3} -p 80 -t {threads} -o {output.out3}
"""

# STATS #

rule fastk:
    input:
        reads=config["reads"],
        hic1=config["hic1"],
        hic2=config["hic2"]
    output:
        ktab=touch("RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".ktab"),
        hist=touch("RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".hist")
    threads:
        workflow.cores
    priority: 1
    params:
        out="RESULTS/STATISTICS_POL/FASTK/" + PREFIX,
        kmers=config["modules"]["fastk"]["kmer"]
    envmodules:
        config["modules"]["fastk"]["path"]
    shell: """
        FastK -v -t1 -T{threads} -k{params.kmers} -N{params.out} {input.reads} {input.hic1} {input.hic2}
    """

rule histex:
    input:
        "RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".hist"
    output:
        "RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".hist.txt"
    priority: 1
    envmodules:
        config["modules"]["fastk"]["path"]
    shell: """
        Histex -G {input} > {output}
    """

rule MerquryFK:
    input:
        hist="RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".hist",
        ktab="RESULTS/STATISTICS_POL/FASTK/" + PREFIX + ".ktab",
        asm=lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/FASTK/{basename}.completeness.stats"
    threads:
        workflow.cores
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        merqury_dir=config["modules"]["merqury"]["path"],
        my_prefix=PREFIX,
        outprefix="{basename}"
    envmodules:
        config["modules"]["fastk"]["path"]
    shell: """
        export PATH=$PATH:{params.merqury_dir} 
        export PATH=$PATH:{params.merqury_dir}/*
        cd RESULTS/STATISTICS_POL/FASTK/
        MerquryFK -T{threads} {params.my_prefix} ../../../{input.asm} {params.outprefix}
    """


rule seqkit_assembly_stats:
    input:
        lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/SEQKIT/{basename}.assemblyStat.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
    """

rule gfastats:
    input:
        fasta=lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/GFASTATS/{basename}_gfastats_report.txt"
    conda:
        "../envs/gfastas.yaml"
    threads:
        min(workflow.cores,10)
    priority: 1
    params:
        genome_size=config["genome_size"]
    shell: """
        gfastats \
        {input.fasta} \
        {params.genome_size} \
        --tabular \
        --threads {threads} \
        --nstar-report \
        --seq-report \
        --out-size scaffolds > {output}
    """

rule quast:
    input:
        contigs=lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
    output:
        "RESULTS/STATISTICS_POL/QUAST/{basename}/quast.log"
    threads: 
        min(workflow.cores,20)
    priority: 1
    params:
        "RESULTS/STATISTICS_POL/QUAST/{basename}"
    conda:
        "../envs/quast.yaml"
    shell: """
        quast.py \
        --eukaryote \
        --large \
        -t {threads} \
        --min-contig 300 \
        -o {params} \
        {input.contigs}
    """

rule busco:
    input:
        lambda wildcards: get_input(wildcards, FA_FILES, BASENAMES, 'basename')
    output:
        summary="RESULTS/BUSCO/{basename}/short_summary.specific.{db}.{basename}.txt",
        missing="RESULTS/BUSCO/{basename}/run_{db}/missing_busco_list.tsv"
    threads:
        min(workflow.cores,20)
    priority: 1
    params:
        dataset_dir=config["busco"]["dataset_folder"],
        out_dir="RESULTS/BUSCO/",
        run_name="{basename}",
        lineage=config["busco"]["lineage"],
        plot_wd="RESULTS/BUSCO/{basename}/"
    conda:
        "../envs/busco.yaml"
    shell: """
        busco -m genome -f \
        -c {threads} \
        -i {input} \
        -o {params.run_name} \
        --download_path {params.dataset_dir} \
        --out_path {params.out_dir} \
        -l {params.lineage} \
        --offline
    """

 
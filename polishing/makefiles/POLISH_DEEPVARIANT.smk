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
    # polishing
        "RESULTS/ASM_POLISHING_DEEPVARIANT/FINAL/" + run + ".fasta",sample=samples,
    # stats
        # expand("RESULTS/BUSCO/{sample}_" + run + "/run_{db}/missing_busco_list.tsv", db = config["busco"]["lineage"], sample = samples),
        # expand("RESULTS/BUSCO/{sample}_" + run + "/short_summary.specific.{db}.{sample}_" + run + ".txt", db = config["busco"]["lineage"], sample = samples),
        # expand("RESULTS/STATISTICS_POL/GFASTATS/{sample}_" + run + "_gfastats_report.txt", sample = samples),
        # expand("RESULTS/STATISTICS_POL/QUAST/{sample}_" + run + "/quast.log", sample = samples),
        # expand("RESULTS/STATISTICS_POL/FASTK/{sample}_" + run + ".completeness.stats", sample = samples),
        # expand("RESULTS/STATISTICS_POL/SEQKIT/{sample}_" + run + ".assemblyStat.txt", sample = samples)



# Function:

def get_chromosomes():
	with open('RESULTS/ASM_POLISHING_DEEPVARIANT/" + run + "/interval.list', 'r') as file:
		return file.read().splitlines()


chromosomes = get_chromosomes()


#replace commas and spaces with underscores
rule correct_headers:
    input:
        lambda wildcards: next(os.path.join(asm_path, f"{wildcards.sample}.{ext}")
                               for ext in extensions
                               if glob.glob(os.path.join(asm_path, f"{wildcards.sample}.{ext}")))
    output:
        touch("RESULTS/CORRECTED_HEADERS/" + run + "_corrected_headers.fasta")
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit replace -p '[ ,]' -r '_' {input} > {output}
    """

rule index_ref:
    input:
        "RESULTS/CORRECTED_HEADERS/" + run + "_corrected_headers.fasta"
    output:
        fai="RESULTS/CORRECTED_HEADERS/" + run + "_corrected_headers.fasta.fai"
    conda:
        "../envs/samtools_minimap.yaml"
    shell: """
        samtools faidx {input}
    """

rule align_reads:
    input:
        ref=lambda wildcards: next(os.path.join(asm_path, f"{wildcards.sample}.{ext}")
                               for ext in extensions
                               if glob.glob(os.path.join(asm_path, f"{wildcards.sample}.{ext}")))
    output:
        bam="RESULTS/ASM_POLISHING_DEEPVARIANT/alignments/" + run + ".bam"
    threads:
        min(workflow.cores,20)
    conda:
        "../envs/pbmm2.yaml"
    params:
        prefix=PREFIX,
        pbmm2_ops=config["pbmm2"]["ops"],
        reads=config["reads"]
    shell: """
        pbmm2 align {input.ref} \
        {params.reads} {output.bam} \
        {params.pbmm2_ops} \
        -j {threads} \
        --rg '@RG\\tID:{params.prefix}\\tSM:{params.prefix}'
    """

# split genome into chromosomes
# two-pass: read file twice to reduce memory usage
# split by the full id name of each chromosome
checkpoint split_by_chromosomes:
    input:
        fasta="RESULTS/CORRECTED_HEADERS/" + run + "_corrected_headers.fasta",
        fai="RESULTS/CORRECTED_HEADERS/" + run + "_corrected_headers.fasta.fai"
    output:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.fasta"
    conda:
        "../envs/seqkit.yaml"
    params:
        outdir="RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/"
    shell:
        """
        seqkit split --quiet --two-pass -i --by-id-prefix "" --out-dir {params.outdir} {input.fasta}
        """


# delete if not useful
checkpoint get_chromosomes:
    input: 
        "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.fasta"
    output:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/" + run + "/interval.list"
    params:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/"
    shell: """
        ls {params} > {output}
    """

rule index_chromosomes:
    input:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.fasta"
    output:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.fasta.fai"
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools faidx {input}
    """

rule create_bed:
    input:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.fasta.fai"
    output:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.bed"
    shell: """
        awk '{{printf("%s\\t0\\t%s\\n",$1,$2)}}' {input}  > {output}
    """


rule split_aligned_reads:
    input:
        bam="RESULTS/ASM_POLISHING_DEEPVARIANT/alignments/" + run + ".bam",
        bed="RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.bed"
    output:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/alignments_split/" + run + "_{chrom}_filt.bam"
    threads:
        min(workflow.cores,20)
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools view -@{threads} -bh -F 2308 -M -L {input.bed} {input.bam} > {output}
        samtools index {output}
    """

# run dv per sample per chromosome
rule run_deepvariant:
    input:
        ref = "RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.fasta",
        bam = "RESULTS/ASM_POLISHING_DEEPVARIANT/alignments_split/" + run + "_{chrom}_filt.bam"
    output:
        vcf = "RESULTS/ASM_POLISHING_DEEPVARIANT/deepvariant/" + run + "_{chrom}/all.filt.dv.vcf"
    threads: 
        min(workflow.cores,20)
    container:
        "/share/scientific_bin/singularity/containers/deepvariant_1.2.0.sif"
    params:
        "output/intermediate_files/"
    shell: """ 
        run_deepvariant \
        --model_type=PACBIO \
        --ref={input.ref} \
        --reads={input.bam} \
        --output_vcf={output.vcf} \
        --intermediate_results_dir {params} \
        --num_shards={threads}
    """

rule filter_deepvariant:
    input:
        chrom="RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.fasta",
        vcf="RESULTS/ASM_POLISHING_DEEPVARIANT/deepvariant/" + run + "_{chrom}/all.filt.dv.vcf"
    output:
        vcf="RESULTS/ASM_POLISHING_DEEPVARIANT/deepvariant_filtered/" + run + "_{chrom}/all.filt.dv.filt.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell: """
        bcftools view -f 'PASS' -i 'GT="1/1"' {input.vcf} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
    """

rule create_consensus_fasta:
    input:
        fasta="RESULTS/ASM_POLISHING_DEEPVARIANT/genome_split/chromosomes/" + run + "/{chrom}.fasta",
        vcf="RESULTS/ASM_POLISHING_DEEPVARIANT/deepvariant_filtered/" + run + "_{chrom}/all.filt.dv.filt.vcf.gz"
    output:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/chrom_consensus/" + run + "_{chrom}.consensus.fasta"
    conda:
        "../envs/bcftools.yaml"
    shell: """
        bcftools consensus -f {input.fasta} {input.vcf} > {output}
    """

# Assuming you have a list or way to dynamically fetch chromosomes after the checkpoint
# will snakemake understand that this is supposed to happen per sample only?

# Function to get paths for all chromosome-specific consensus FASTA files for merging
def get_chromosomes_for_sample(wildcards):
    interval_list_path = f"RESULTS/ASM_POLISHING_DEEPVARIANT/{wildcards.sample}_{run}/interval.list"
    try:
        with open(interval_list_path, 'r') as file:
            chroms = file.read().splitlines()
        return chroms
    except FileNotFoundError:
        return []

# Rule to merge chromosome-specific consensus FASTA files into a single FASTA file per sample
rule merge_chroms:
    input: 
        get_chromosomes_for_sample
    output: 
        expand("RESULTS/ASM_POLISHING_DEEPVARIANT/FINAL/" + run + ".fasta")
    shell: 
        "cat {input} > {output}"




# rule merge_chroms:
#     input:
#         "RESULTS/ASM_POLISHING_DEEPVARIANT/chrom_consensus/{sample}_" + run + "_{chrom}.consensus.fasta"
#     output:
#         "RESULTS/ASM_POLISHING_DEEPVARIANT/FINAL/{sample}_" + run + ".fasta"
#     shell:
#         "cat {input} > {output}"


# STATS

rule fastk:
    input:
        config["reads"]
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
        FastK -v -t1 -T{threads} -k{params.kmers} -N{params.out} {input}
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
        asm="RESULTS/ASM_POLISHING_DEEPVARIANT/FINAL/{sample}_" + run + ".fasta"
    output:
        "RESULTS/STATISTICS_POL/FASTK/{sample}_" + run + ".completeness.stats"
    threads:
        workflow.cores
    priority: 1
    conda:
        "../envs/merqury.yaml"
    params:
        merqury_dir=config["modules"]["merqury"]["path"],
        my_prefix=PREFIX
    envmodules:
        config["modules"]["fastk"]["path"]
    shell: """
        export PATH=$PATH:{params.merqury_dir} 
        export PATH=$PATH:{params.merqury_dir}/*
        cd RESULTS/STATISTICS_POL/FASTK/
        MerquryFK -T{threads} {params.my_prefix} ../../../{input.asm} {wildcards.sample}
    """


rule seqkit_assembly_stats:
    input:
        "RESULTS/ASM_POLISHING_DEEPVARIANT/FINAL/{sample}_" + run + ".fasta"
    output:
        "RESULTS/STATISTICS_POL/SEQKIT/{sample}_" + run + ".assemblyStat.txt"
    priority: 1
    conda:
        "../envs/seqkit.yaml"
    shell: """
        seqkit stats -a {input} > {output}
    """

rule gfastats:
    input:
        fasta="RESULTS/ASM_POLISHING_DEEPVARIANT/FINAL/{sample}_" + run + ".fasta"
    output:
        "RESULTS/STATISTICS_POL/GFASTATS/{sample}_" + run + "_gfastats_report.txt"
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
        contigs = "RESULTS/ASM_POLISHING_DEEPVARIANT/FINAL/{sample}_" + run + ".fasta"
    output:
        "RESULTS/STATISTICS_POL/QUAST/{sample}_" + run + "/quast.log"
    threads: 
        min(workflow.cores,20)
    priority: 1
    params:
        "RESULTS/STATISTICS_POL/QUAST/{sample}_" + run
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
        "RESULTS/ASM_POLISHING_DEEPVARIANT/FINAL/{sample}_" + run + ".fasta"
    output:
        summary="RESULTS/BUSCO/{sample}_" + run + "/short_summary.specific.{db}.{sample}_" + run + ".txt",
        missing="RESULTS/BUSCO/{sample}_" + run + "/run_{db}/missing_busco_list.tsv"
    threads:
        min(workflow.cores,20)
    priority: 1
    params:
        dataset_dir=config["busco"]["dataset_folder"],
        out_dir="RESULTS/BUSCO/",
        run_name="{sample}_" + run,
        lineage=config["busco"]["lineage"],
        plot_wd="RESULTS/BUSCO/{sample}_" + run + "/"
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

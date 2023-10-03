############ Purge erroneous duplications ############


# Programs and Descriptions in current file:
# Split the FASTA file wherever 'N' nucleotides occur.
#   - Program: Custom script (`split_fa`) from purge_dups package.

# - Map the assembly to itself to identify regions of similarity.
#   - Program: Minimap2.

# - Map the split genome to HiFi reads.
#   - Program: Minimap2.

# - Calculate the base-level coverage for the mapped reads.
#   - Program: `pbcstat`.

# - Determine coverage thresholds to differentiate between haploid and diploid regions.
#   - Program: `calcuts`.

# - Purge Duplicates
#   - Program: `purge_dups`.

# - Extract sequences that are not marked as duplicates with get_seqs.
#   - Program: `get_seqs` from purgedups package.

# - Index the purged sequences
#   - Program: Samtools (`samtools faidx`).


# split fasta file where Ns occur

rule SPLIT_AT_Ns:
    input:
        lambda wildcards: get_input(wildcards, HIFIASM_FILES, HIFIASM_BASENAMES, 'basename')
    output:
        temp("RESULTS/PURGE_DUPS/{basename}_split.fasta")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{basename}.split_fa.log"
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        split_fa {input} > {output}
    """

# xasm5=asm-to-ref mapping, for ~0.1 sequence divergence (assembly to self)
# -I=split index for every ~NUM input bases

rule MINIMAP2_GENOME_TO_SELF:
    input:
        "RESULTS/PURGE_DUPS/{basename}_split.fasta"
    output:
        temp("RESULTS/MINIMAP2/{basename}_genome_split.paf")
    threads:
        workflow.cores
    params:
        i_split=config["INDEX_SPLIT"]
    conda:
        "../envs/minimap2.yaml"
    shell: """
        minimap2 -I {params.i_split} -t {threads} -xasm5 -DP {input} {input} > {output}
    """

rule MINIMAP2_GENOME_TO_READS:
    input:
        splitf="RESULTS/PURGE_DUPS/{basename}_split.fasta",
        hifi_reads=lambda wildcards: "DATA/" + PREFIX + ".fastq.gz" if SUBREADS else config["HIFI"]
    output:
        temp("RESULTS/MINIMAP2/{basename}_reads_split.paf")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{basename}.minimap2b.log"
    params:
        i_split=config["INDEX_SPLIT"]
    conda:
        "../envs/minimap2.yaml"
    shell: """
        minimap2 -I {params.i_split} -t {threads} -x map-pb {input.splitf} {input.hifi_reads} > {output}
    """

#Calculate haploid/diploid coverage threshold and remove haplotype duplicates from assembly

rule PBCSTAT:
    input:
        "RESULTS/MINIMAP2/{basename}_reads_split.paf"
    output:
        basecov="RESULTS/PURGE_DUPS/{basename}/PB.base.cov",
        stat="RESULTS/PURGE_DUPS/{basename}/PB.stat" #, PB.cov.wig
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{basename}.pbcstat.log"
    params:
        prefix="RESULTS/PURGE_DUPS/{basename}"
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        pbcstat -O {params.prefix} {input}
    """

rule CALCUTS:
    input:
        "RESULTS/PURGE_DUPS/{basename}/PB.stat"
    output:
        "RESULTS/PURGE_DUPS/{basename}/cutoffs"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{basename}.calcuts.log"
    params:
        prefix="RESULTS/PURGE_DUPS/{basename}"
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        calcuts {input} > {output}
    """

rule PURGE_DUPS:
    input:
        basecov="RESULTS/PURGE_DUPS/{basename}/PB.base.cov",
        cutoffs="RESULTS/PURGE_DUPS/{basename}/cutoffs",
        genome_paf="RESULTS/MINIMAP2/{basename}_genome_split.paf"
    output:
        temp("RESULTS/PURGE_DUPS/{basename}_dups.bed")
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{basename}.purge_dups.log"
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        purge_dups -2 -c {input.basecov} -T {input.cutoffs} {input.genome_paf} > {output}
    """

rule GET_SEQS:
    input:
        purge_duped="RESULTS/PURGE_DUPS/{basename}_dups.bed",
        fasta="RESULTS/GENOME_ASSEMBLY/{basename}.fasta"
    output:
        "RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa"
    threads:
        workflow.cores
    log:
        "RESULTS/LOG/{basename}.get_seqs.log"
    params:
        prefix="RESULTS/PURGE_DUPS/{basename}_seqs_purged"
    conda:
        "../envs/purge_dups.yaml"
    shell: """
        get_seqs -e -p {params.prefix} {input.purge_duped} {input.fasta}
    """


rule INDEX_SEQS:
    input:
        "RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa"
    output:
        "RESULTS/PURGE_DUPS/{basename}_seqs_purged.hap.fa.fai"
    conda:
        "../envs/samtools.yaml"
    shell: """
        samtools faidx {input}
    """
